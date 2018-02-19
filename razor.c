#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <getopt.h>
#include <stdbool.h>
#include <zlib.h>

#define VERSION 0.1
#define MAX_SEQ_SIZE 32786
#define MAX_HEADER_SIZE 4096
#define MAX_ADAPTERS 64
#define MAX_ADAPTER_LEN 256

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

struct Record {
	char* header;
	char* sequence;
	char* comment;
	char* quality;
	int len;
};

struct TrimArgs {
	char* adapter;
	char* adapter_file;
	int quality_threshold;
	int min_length;
	int min_mean_quality;
	int min_adapter_match;
	int phred_encoding;
	bool keep_empty;
};

struct Adapters {
	char* seqs[MAX_ADAPTERS];
	int* fail_mtrx[MAX_ADAPTERS];
	int lengths[MAX_ADAPTERS];
};

struct StatsReport {
	unsigned long reads_read;
	unsigned long reads_quality_filtered;
	unsigned long reads_adapter_filtered;
	unsigned long reads_adapter_trimmed;
	unsigned long reads_quality_trimmed;
	unsigned long bases_adapter_trimmed;
	unsigned long bases_quality_trimmed;
};


static inline void* zmalloc(size_t size)
{
	void* ptr = malloc(size);
	if(!ptr && size)
	{
		fprintf(stderr, "Could not allocate %zu bytes\n", size);
		exit(EXIT_FAILURE);
	}
	return ptr;
}



static inline void* zcalloc(size_t nmemb, size_t size)
{
	void* ptr = calloc(nmemb, size);
	if(!ptr && nmemb && size)
	{
		fprintf(stderr, "Could not allocate %zu bytes\n", nmemb * size);
		exit(EXIT_FAILURE);
	}
	return ptr;
}



int num_places(unsigned long n) {
	if (n < 10) return 1;
	return 1 + num_places (n / 10);
}



void print_stats(struct StatsReport stats, bool trim_adapters, bool trim_quality)
{
	unsigned long max;
	max = MAX(stats.bases_quality_trimmed, stats.bases_adapter_trimmed);
	max = MAX(max, stats.reads_read);
	int max_width = num_places(max);

	fprintf(stderr, "%*lu reads were read\n", max_width, stats.reads_read);
	if (trim_adapters){
		fprintf(stderr, "%*lu reads underwent adapter trimming\n", max_width, stats.reads_adapter_trimmed);
		fprintf(stderr, "%*lu of these were filtered\n", max_width, stats.reads_adapter_filtered);
		fprintf(stderr, "%*lu adapter bases were trimmed\n", max_width, stats.bases_adapter_trimmed);
	}

	if (trim_quality) {
		fprintf(stderr, "%*lu reads underwent quality trimming\n", max_width, stats.reads_quality_trimmed);
		fprintf(stderr, "%*lu of these were filtered\n", max_width, stats.reads_quality_filtered);
		fprintf(stderr, "%*lu bases were quality trimmed\n", max_width, stats.bases_quality_trimmed);
	}
}



static inline struct Record* init_record()
{
	struct Record* rec = zmalloc(sizeof(struct Record));
	rec->header = zmalloc(2 * MAX_HEADER_SIZE + 2 * MAX_SEQ_SIZE);
	rec->sequence = rec->header + MAX_HEADER_SIZE;
	rec->comment = rec->sequence + MAX_SEQ_SIZE;
	rec->quality = rec->comment + MAX_HEADER_SIZE;
	return rec;
}



bool quality_trim(int qual, int phred, int min_length, int min_mean_qual, 
		   struct Record *rec, struct StatsReport* stats)
{
	int start = 0, end = 0, sum = 0, max = 0, i = 0;
	int q = qual + phred;

	while (rec->quality[i] < q && i < rec->len) i++;
	start = i;
	while (i < rec->len) {
		while (rec->quality[i] >= q && i < rec->len) {
			sum += rec->quality[i++] - q;
		}
		if (sum > max) {
			max = sum;
			end = i;
		}
		while (rec->quality[i] < q && i < rec->len) {
			sum += rec->quality[i++] - q;
		}
	}

	rec->quality += start;
	rec->sequence += start;
	rec->quality[end]  = '\0';
	rec->sequence[end] = '\0';
	int length = end - start;

	if (rec->len == length) {
		return true;
	} else {
		stats->reads_quality_trimmed++;
		stats->bases_quality_trimmed += rec->len - length;
		if ((float)max / (float)length < min_mean_qual ||
			length < min_length) {
			rec->sequence[0] = '\0';
			rec->quality[0] = '\0';
			stats->reads_quality_filtered++;
			return false;
		}
		return true;
	}

}



bool adapter_trim(int min_length, int num_adapters, struct Adapters* adapters,
				  struct Record* rec, struct StatsReport* stats)
{
	int i, j;
	for (int k = 0; k < num_adapters; k++) {
		i = 0, j = 0;
		while (i < rec->len) {
			if (rec->sequence[i] == adapters->seqs[k][j]) {
				if (j == adapters->lengths[k] - 1) {
					stats->reads_adapter_trimmed++;
					stats->bases_adapter_trimmed += rec->len - i;
					rec->len = i - j;
					if (i - j < min_length) {
						rec->sequence[0] = '\0';
						rec->quality[0]  = '\0';
						stats->reads_adapter_filtered++;
						return false;
					} else {
						rec->sequence[i - j] = '\0';
						rec->quality[i - j]  = '\0';
						return true;
					}
				} else {
					i++;
					j++;
				}
			} else if (j > 0) {
				j = adapters->fail_mtrx[k][j - 1];
			} else {
				i++;
			}
		}

	}
	return true;
}



void kmp_preprocess(struct Adapters* adapters, int num_adapters)
{
	int i, j;
	for (int k = 0; k < num_adapters; k++){
		i = 1, j = 0;
		while (i < adapters->lengths[k]){
			if (adapters->seqs[k][i] == adapters->seqs[k][j]){
				adapters->fail_mtrx[k][i] = j + 1;
				i++; j++;
			} else if (j > 0)
				j = adapters->fail_mtrx[k][j-1];
			else {
				adapters->fail_mtrx[k][i] = 0;
				i++;
			}
		}
	}
}



static inline int parse_fasta(gzFile input, int min_match, char* seq)
{
	if (!gzgets(input, seq, MAX_ADAPTER_LEN) || gzeof(input))
		return 0;
	if (seq[0] != '>') {
		fprintf(stderr, "Expected fasta header starting with '>'; instead got '%c'.\n", seq[0]);
		exit(EXIT_FAILURE);
	}

	int i = 0;
	gzgets(input, seq, MAX_ADAPTER_LEN);
	while (seq[i] != '\n'){
		switch(seq[i]){
			case 'A':
			case 'C':
			case 'T':
			case 'G':
			case 'U':
			case 'N':
				break;
			default:
				fprintf(stderr, "Invalid character '%c' detected in adapter file\n"
								"Valid characters: A,C,G,T,U,N\n",
						seq[i]);
				exit(EXIT_FAILURE);
		}
		i++;
	}

	int length = MIN(i, min_match);
	seq[length] = '\0';
	return length;
}



int init_adapters(struct Adapters* adapters, char* adapter_file, char* adapter, int min_match)
{
	adapters->seqs[0] = zcalloc(MAX_ADAPTERS * MAX_ADAPTER_LEN, 1);
	adapters->fail_mtrx[0] = zcalloc(MAX_ADAPTERS * MAX_ADAPTER_LEN, sizeof(int));

	int i = 0;

	if (adapter){
		adapter[MIN(min_match, strlen(adapter))] = '\0';
		adapters->lengths[i] = strlen(adapter);
		strcpy(adapters->seqs[i], adapter);
		i++;
	}

	if (adapter_file) {
		gzFile input = gzopen(adapter_file, "r");
		if (!input) {
			fprintf(stderr, "Error opening '%s': %s\n", adapter_file, strerror(errno));
			exit(EXIT_FAILURE);
		}

		while (i < MAX_ADAPTERS) {
			if ((adapters->lengths[i] = parse_fasta(input, min_match, adapters->seqs[i])) <= 0)
				break;
			adapters->seqs[i + 1] = adapters->seqs[i] + adapters->lengths[i];
			adapters->fail_mtrx[i + 1] = adapters->fail_mtrx[i] + adapters->lengths[i];
			i++;
		}
		gzclose_r(input);
	}

	adapters->seqs[0] = realloc(adapters->seqs[0], i * MAX_ADAPTER_LEN);
	adapters->fail_mtrx[0] = realloc(adapters->fail_mtrx[0], i * MAX_ADAPTER_LEN);
	return i;
}



static inline bool gzgetrec(gzFile input, struct Record* rec)
{
	if (!gzgets(input, rec->header, MAX_HEADER_SIZE)) return 0;
	gzgets(input, rec->sequence, MAX_SEQ_SIZE);
	gzgets(input, rec->comment, MAX_HEADER_SIZE);
	gzgets(input, rec->quality, MAX_SEQ_SIZE);

	rec->len = strlen(rec->sequence) - 1;

	if (rec->len != strlen(rec->quality) - 1){
		fprintf(stderr, "Error: lengths of sequence and quality strings "
						"unequal:\n%s%s%s%s\n",
						rec->header, rec->sequence,
						rec->comment, rec->quality);
		exit(EXIT_FAILURE);
	}

	if (rec->header[0] != '@'){
		fprintf(stderr, "Malformed fastq: header should start with '@', not"
			"'%c':\n%s\n", rec->header[0], rec->header);
		exit(EXIT_FAILURE);
	}
	if (rec->comment[0] != '+'){
		fprintf(stderr, "Malformed fastq: comment should start with '+', not"
			"'%c':\n%s\n", rec->comment[0], rec->comment);
		exit(EXIT_FAILURE);
	}

	rec->sequence[rec->len] = '\0';
	rec->quality[rec->len] = '\0';
	return 1;
}



void process_reads(const char* infile, const char* outfile,
				   bool trim_adapters, bool trim_quality,
				   const struct TrimArgs* trim_args)
{
	gzFile input;
	if (!strcmp(infile, "-"))
		input = gzdopen(STDIN_FILENO, "r");
	else
		input = gzopen(infile, "r");
	if (!input){
		fprintf(stderr, "Error opening '%s': %s\n", infile, strerror(errno));
		exit(EXIT_FAILURE);
	}

	bool gzip_output = false;
	void* output;
	if (!outfile || !strcmp(outfile, "-")) {
		output = fdopen(STDOUT_FILENO, "w");
	}
	else {
		char* dot = strrchr(outfile, '.');
		if (dot && !strcmp(dot, ".gz")) {
			output = gzopen(outfile, "w");
			gzip_output = true;
		}
		else {
			output = fopen(outfile, "w");
		}
	}
	if (!output){
		fprintf(stderr, "Error opening '%s': %s\n", outfile, strerror(errno));
		exit(EXIT_FAILURE);
	}

	struct Adapters adapters;
	int num_adapters = 0;
	if (trim_adapters){
		num_adapters = init_adapters(&adapters, trim_args->adapter_file,
									 			trim_args->adapter,
									 			trim_args->min_adapter_match);

		kmp_preprocess(&adapters, num_adapters);
	}

	struct Record* rec_ptr = init_record();
	struct Record rec;
	struct StatsReport stats = {0};
	bool passed;
	while (true) {
		if (gzeof(input)) break;
		if (!gzgetrec(input, rec_ptr)) break;
		stats.reads_read++;

		rec = *rec_ptr;
		passed  = true;

		if (trim_adapters)
			passed = adapter_trim(trim_args->min_length, num_adapters, &adapters,
								  &rec, &stats);

		if (passed && trim_quality)
			passed = quality_trim(trim_args->quality_threshold,
								  trim_args->phred_encoding,
								  trim_args->min_length,
								  trim_args->min_mean_quality,
								  &rec, &stats);

		if ( passed || trim_args->keep_empty) {
			if (gzip_output) {
				if (!gzprintf(output, "%s%s\n%s%s\n", rec.header, rec.sequence,
							  rec.comment, rec.quality)) {
					fprintf(stderr, "Error writing to '%s': %s\n", outfile,
							strerror(errno));
					exit(EXIT_FAILURE);
				}
			} else {
				if (!fprintf(output, "%s%s\n%s%s\n", rec.header, rec.sequence,
							 rec.comment, rec.quality)) {
					fprintf(stderr, "Error writing to '%s': %s\n", outfile,
							strerror(errno));
					exit(EXIT_FAILURE);
				}
			}

		}

	}
	print_stats(stats, trim_adapters, trim_quality);

	gzclose_r(input);
	if (gzip_output)
		gzclose_w(output);
	else
		fclose(output);
//	free(rec_ptr->header);
//	free(rec_ptr);
//	free(adapters.seqs[0]);
//	free(adapters.fail_mtrx[0]);
}



void print_usage()
{
	fprintf(stderr, "Usage: razor [options] [in.fq]\n");
	fprintf(stderr, "Options: \n");
	fprintf(stderr, "  -o, --output=FILE		The file to output processed reads to. If unspecified, reads are written to stdout.\n");
	fprintf(stderr, "  -q, --quality-threshold=INT	Trim bases from 5' and 3' ends using the ERNE-FILTER algorithm.\n");
	fprintf(stderr, "  -l, --min-length=INT		After processing, discard sequences below a minimum length [20].\n");
	fprintf(stderr, "  -m, --min-mean-quality=INT	After quality trimming, discard reads below a minimum mean quality.\n");
	fprintf(stderr, "  -a, --adapter=STR		Adapter to be trimmed from the 3' ends of the reads.\n");
	fprintf(stderr, "  -f, --adapter-file=FILE	Fasta file containing adapter sequences.\n");
	fprintf(stderr, "  -M, --min-adapter-match=INT	The minimum number of leading bases a subsequence must share with an adapter\n"
	    		"				for a match to occur [12].\n");
	fprintf(stderr, "  -k, --keep-empty		Sequences filtered according to the min-length and min-mean-quality options\n"
			"				have their sequence and quality strings truncated to empty and are retained.\n");
	fprintf(stderr, "  -p, --phred=[33|64]		The phred encoding of the quality scores in the input fastq [33].\n");
	fprintf(stderr, "      --version			Display version information and exit.\n");
	fprintf(stderr, "      --help			Display this help text and exit.\n");

}



void print_version()
{
	fprintf(stderr, "Razor version %.1f\n", VERSION);
}



static struct option long_opts[] = {
	{"output", required_argument, 0, 'o'},
	{"quality-threshold", required_argument, 0, 'q'},
	{"min-length", required_argument, 0, 'l'},
	{"min-mean-quality", required_argument, NULL, 'm'},
	{"adapter", required_argument, 0, 'a'},
	{"adapter-file", required_argument, 0, 'f'},
	{"min-adapter-match", required_argument, 0, 'M'},
	{"keep-empty", no_argument, 0, 'k'},
	{"phred", required_argument, 0, 'p'},
	{"help", no_argument, NULL, 2},
	{"version", no_argument, NULL, 3},
	{0, 0, 0, 0}
};



int main(int argc, char *argv[])
{
	static char* infile;
	static char* outfile;
	static bool trim_adapters = false;
	static bool trim_quality = false;
	static struct TrimArgs trim_args = {
		.quality_threshold = 0,
		.min_length = 0,
		.min_mean_quality = 0,
		.adapter = NULL,
		.adapter_file = NULL,
		.min_adapter_match = 0,
		.keep_empty = false,
		.phred_encoding = 33,
	};

	static const char* short_opts = "o:q:l:m:a:f:M:p:k";
	int c;
	int index = 0;
	while((c = getopt_long(argc, argv, short_opts, long_opts, &index)) != -1){
		switch(c){
			case 'o': outfile = optarg; break;
			case 'q': trim_args.quality_threshold = atoi(optarg); break;
			case 'l': trim_args.min_length = atoi(optarg); break;
			case 'm': trim_args.min_mean_quality = atoi(optarg); break;
			case 'a': trim_args.adapter = optarg; break;
			case 'f': trim_args.adapter_file = optarg; break;
			case 'M': trim_args.min_adapter_match = atoi(optarg); break;
			case 'k': trim_args.keep_empty = true; break;
			case 'p':
				if (atoi(optarg) == 33 || atoi(optarg) == 64){
					trim_args.phred_encoding = atoi(optarg);
					break;
				}
				print_usage();
				return 1;
			case 2: print_usage(); return 0;
			case 3: print_version(); return 0;
			default:
				print_usage();
				return 1;
		}
	}

	if (argc - optind != 1) {
		print_usage();
		return EXIT_FAILURE;
	}

	if (trim_args.adapter || trim_args.adapter_file) {
		trim_adapters = true;
		if (!trim_args.min_adapter_match)
			trim_args.min_adapter_match = 12;
	} else {
		if (trim_args.min_adapter_match){
			fprintf(stderr, "Error: option --min-adapter-match specified without adapter sequence(s)\n");
			return EXIT_FAILURE;
		}
//		trim_adapters = false;
	}

	if (trim_args.quality_threshold)
		trim_quality = true;
	else {
		if (trim_args.min_mean_quality){
			fprintf(stderr, "Error: option --min-mean-quality requires the -q option\n");
			return EXIT_FAILURE;
		}
		if (!trim_adapters){
			fprintf(stderr, "Error: must specify at least one of the following options: -q, -a, -f\n");
			return EXIT_FAILURE;
		}
//		trim_quality = false;
	}

	if (!trim_args.min_length)
		trim_args.min_length = 20;

	infile = argv[optind];
	process_reads(infile, outfile, trim_adapters, trim_quality, &trim_args);

	return EXIT_SUCCESS;
}


