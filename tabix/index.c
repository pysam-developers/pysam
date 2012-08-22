#include <ctype.h>
#include <assert.h>
#include <sys/stat.h>
#include "khash.h"
#include "ksort.h"
#include "kstring.h"
#include "bam_endian.h"
#ifdef _USE_KNETFILE
#include "knetfile.h"
#endif
#include "tabix.h"

#define TAD_MIN_CHUNK_GAP 32768
// 1<<14 is the size of minimum bin.
#define TAD_LIDX_SHIFT    14

typedef struct {
	uint64_t u, v;
} pair64_t;

#define pair64_lt(a,b) ((a).u < (b).u)
KSORT_INIT(offt, pair64_t, pair64_lt)

typedef struct {
	uint32_t m, n;
	pair64_t *list;
} ti_binlist_t;

typedef struct {
	int32_t n, m;
	uint64_t *offset;
} ti_lidx_t;

KHASH_MAP_INIT_INT(i, ti_binlist_t)
KHASH_MAP_INIT_STR(s, int)

struct __ti_index_t {
	ti_conf_t conf;
	int32_t n, max;
	khash_t(s) *tname;
	khash_t(i) **index;
	ti_lidx_t *index2;
};

struct __ti_iter_t {
	int from_first; // read from the first record; no random access
	int tid, beg, end, n_off, i, finished;
	uint64_t curr_off;
	kstring_t str;
	const ti_index_t *idx;
	pair64_t *off;
};

typedef struct {
	int tid, beg, end, bin;
} ti_intv_t;

ti_conf_t ti_conf_gff = { 0, 1, 4, 5, '#', 0 };
ti_conf_t ti_conf_bed = { TI_FLAG_UCSC, 1, 2, 3, '#', 0 };
ti_conf_t ti_conf_psltbl = { TI_FLAG_UCSC, 15, 17, 18, '#', 0 };
ti_conf_t ti_conf_sam = { TI_PRESET_SAM, 3, 4, 0, '@', 0 };
ti_conf_t ti_conf_vcf = { TI_PRESET_VCF, 1, 2, 0, '#', 0 };

/***************
 * read a line *
 ***************/

/*
int ti_readline(BGZF *fp, kstring_t *str)
{
	int c, l = 0;
	str->l = 0;
	while ((c = bgzf_getc(fp)) >= 0 && c != '\n') {
		++l;
		if (c != '\r') kputc(c, str);
	}
	if (c < 0 && l == 0) return -1; // end of file
	return str->l;
}
*/

/* Below is a faster implementation largely equivalent to the one
 * commented out above. */
int ti_readline(BGZF *fp, kstring_t *str)
{
	return bgzf_getline(fp, '\n', str);
}

/*************************************
 * get the interval from a data line *
 *************************************/

static inline int ti_reg2bin(uint32_t beg, uint32_t end)
{
	--end;
	if (beg>>14 == end>>14) return 4681 + (beg>>14);
	if (beg>>17 == end>>17) return  585 + (beg>>17);
	if (beg>>20 == end>>20) return   73 + (beg>>20);
	if (beg>>23 == end>>23) return    9 + (beg>>23);
	if (beg>>26 == end>>26) return    1 + (beg>>26);
	return 0;
}

static int get_tid(ti_index_t *idx, const char *ss)
{
	khint_t k;
	int tid;
	k = kh_get(s, idx->tname, ss);
	if (k == kh_end(idx->tname)) { // a new target sequence
		int ret, size;
		// update idx->n, ->max, ->index and ->index2
		if (idx->n == idx->max) {
			idx->max = idx->max? idx->max<<1 : 8;
			idx->index = realloc(idx->index, idx->max * sizeof(void*));
			idx->index2 = realloc(idx->index2, idx->max * sizeof(ti_lidx_t));
		}
		memset(&idx->index2[idx->n], 0, sizeof(ti_lidx_t));
		idx->index[idx->n++] = kh_init(i);
		// update ->tname
		tid = size = kh_size(idx->tname);
		k = kh_put(s, idx->tname, strdup(ss), &ret);
		kh_value(idx->tname, k) = size;
		assert(idx->n == kh_size(idx->tname));
	} else tid = kh_value(idx->tname, k);
	return tid;
}

int ti_get_intv(const ti_conf_t *conf, int len, char *line, ti_interval_t *intv)
{
	int i, b = 0, id = 1, ncols = 0;
	char *s;
	intv->ss = intv->se = 0; intv->beg = intv->end = -1;
	for (i = 0; i <= len; ++i) {
		if (line[i] == '\t' || line[i] == 0) {
            ++ncols;
			if (id == conf->sc) {
				intv->ss = line + b; intv->se = line + i;
			} else if (id == conf->bc) {
				// here ->beg is 0-based.
				intv->beg = intv->end = strtol(line + b, &s, 0);
				if (!(conf->preset&TI_FLAG_UCSC)) --intv->beg;
				else ++intv->end;
				if (intv->beg < 0) intv->beg = 0;
				if (intv->end < 1) intv->end = 1;
			} else {
				if ((conf->preset&0xffff) == TI_PRESET_GENERIC) {
					if (id == conf->ec) intv->end = strtol(line + b, &s, 0);
				} else if ((conf->preset&0xffff) == TI_PRESET_SAM) {
					if (id == 6) { // CIGAR
						int l = 0, op;
						char *t;
						for (s = line + b; s < line + i;) {
							long x = strtol(s, &t, 10);
							op = toupper(*t);
							if (op == 'M' || op == 'D' || op == 'N') l += x;
							s = t + 1;
						}
						if (l == 0) l = 1;
						intv->end = intv->beg + l;
					}
				} else if ((conf->preset&0xffff) == TI_PRESET_VCF) {
					// FIXME: the following is NOT tested and is likely to be buggy
					if (id == 4) {
						if (b < i) intv->end = intv->beg + (i - b);
					} else if (id == 8) { // look for "END="
						int c = line[i];
						line[i] = 0;
						s = strstr(line + b, "END=");
						if (s == line + b) s += 4;
						else if (s) {
							s = strstr(line + b, ";END=");
							if (s) s += 5;
						}
						if (s) intv->end = strtol(s, &s, 0);
						line[i] = c;
					}
				}
			}
			b = i + 1;
			++id;
		}
	}
/*
	if (ncols < conf->sc || ncols < conf->bc || ncols < conf->ec) {
		if (ncols == 1) fprintf(stderr,"[get_intv] Is the file tab-delimited? The line has %d field only: %s\n", ncols, line);
		else fprintf(stderr,"[get_intv] The line has %d field(s) only: %s\n", ncols, line);
		exit(1);
	}
*/
	if (intv->ss == 0 || intv->se == 0 || intv->beg < 0 || intv->end < 0) return -1;
	return 0;
}

static int get_intv(ti_index_t *idx, kstring_t *str, ti_intv_t *intv)
{
	ti_interval_t x;
	intv->tid = intv->beg = intv->end = intv->bin = -1;
	if (ti_get_intv(&idx->conf, str->l, str->s, &x) == 0) {
		int c = *x.se;
		*x.se = '\0'; intv->tid = get_tid(idx, x.ss); *x.se = c;
		intv->beg = x.beg; intv->end = x.end;
		intv->bin = ti_reg2bin(intv->beg, intv->end);
		return (intv->tid >= 0 && intv->beg >= 0 && intv->end >= 0)? 0 : -1;
	} else {
		fprintf(stderr, "[%s] the following line cannot be parsed and skipped: %s\n", __func__, str->s);
		return -1;
	}
}

/************
 * indexing *
 ************/

// requirement: len <= LEN_MASK
static inline void insert_offset(khash_t(i) *h, int bin, uint64_t beg, uint64_t end)
{
	khint_t k;
	ti_binlist_t *l;
	int ret;
	k = kh_put(i, h, bin, &ret);
	l = &kh_value(h, k);
	if (ret) { // not present
		l->m = 1; l->n = 0;
		l->list = (pair64_t*)calloc(l->m, 16);
	}
	if (l->n == l->m) {
		l->m <<= 1;
		l->list = (pair64_t*)realloc(l->list, l->m * 16);
	}
	l->list[l->n].u = beg; l->list[l->n++].v = end;
}

static inline uint64_t insert_offset2(ti_lidx_t *index2, int _beg, int _end, uint64_t offset)
{
	int i, beg, end;
	beg = _beg >> TAD_LIDX_SHIFT;
	end = (_end - 1) >> TAD_LIDX_SHIFT;
	if (index2->m < end + 1) {
		int old_m = index2->m;
		index2->m = end + 1;
		kroundup32(index2->m);
		index2->offset = (uint64_t*)realloc(index2->offset, index2->m * 8);
		memset(index2->offset + old_m, 0, 8 * (index2->m - old_m));
	}
	if (beg == end) {
		if (index2->offset[beg] == 0) index2->offset[beg] = offset;
	} else {
		for (i = beg; i <= end; ++i)
			if (index2->offset[i] == 0) index2->offset[i] = offset;
	}
	if (index2->n < end + 1) index2->n = end + 1;
	return (uint64_t)beg<<32 | end;
}

static void merge_chunks(ti_index_t *idx)
{
	khash_t(i) *index;
	int i, l, m;
	khint_t k;
	for (i = 0; i < idx->n; ++i) {
		index = idx->index[i];
		for (k = kh_begin(index); k != kh_end(index); ++k) {
			ti_binlist_t *p;
			if (!kh_exist(index, k)) continue;
			p = &kh_value(index, k);
			m = 0;
			for (l = 1; l < p->n; ++l) {
				if (p->list[m].v>>16 == p->list[l].u>>16) p->list[m].v = p->list[l].v;
				else p->list[++m] = p->list[l];
			} // ~for(l)
			p->n = m + 1;
		} // ~for(k)
	} // ~for(i)
}

static void fill_missing(ti_index_t *idx)
{
	int i, j;
	for (i = 0; i < idx->n; ++i) {
		ti_lidx_t *idx2 = &idx->index2[i];
		for (j = 1; j < idx2->n; ++j)
			if (idx2->offset[j] == 0)
				idx2->offset[j] = idx2->offset[j-1];
	}
}

ti_index_t *ti_index_core(BGZF *fp, const ti_conf_t *conf)
{
	int ret;
	ti_index_t *idx;
	uint32_t last_bin, save_bin;
	int32_t last_coor, last_tid, save_tid;
	uint64_t save_off, last_off, lineno = 0, offset0 = (uint64_t)-1, tmp;
	kstring_t *str;

	str = calloc(1, sizeof(kstring_t));

	idx = (ti_index_t*)calloc(1, sizeof(ti_index_t));
	idx->conf = *conf;
	idx->n = idx->max = 0;
	idx->tname = kh_init(s);
	idx->index = 0;
	idx->index2 = 0;

	save_bin = save_tid = last_tid = last_bin = 0xffffffffu;
	save_off = last_off = bgzf_tell(fp); last_coor = 0xffffffffu;
	while ((ret = ti_readline(fp, str)) >= 0) {
		ti_intv_t intv;
		++lineno;
		if (lineno <= idx->conf.line_skip || str->s[0] == idx->conf.meta_char) {
			last_off = bgzf_tell(fp);
			continue;
		}
		get_intv(idx, str, &intv);
        if ( intv.beg<0 || intv.end<0 )
        {
            fprintf(stderr,"[ti_index_core] the indexes overlap or are out of bounds\n");
            exit(1);
        }
		if (last_tid != intv.tid) { // change of chromosomes
            if (last_tid>intv.tid )
            {
                fprintf(stderr,"[ti_index_core] the chromosome blocks not continuous at line %llu, is the file sorted? [pos %d]\n",(unsigned long long)lineno,intv.beg+1);
                exit(1);
            }
			last_tid = intv.tid;
			last_bin = 0xffffffffu;
		} else if (last_coor > intv.beg) {
			fprintf(stderr, "[ti_index_core] the file out of order at line %llu\n", (unsigned long long)lineno);
			exit(1);
		}
		tmp = insert_offset2(&idx->index2[intv.tid], intv.beg, intv.end, last_off);
		if (last_off == 0) offset0 = tmp;
		if (intv.bin != last_bin) { // then possibly write the binning index
			if (save_bin != 0xffffffffu) // save_bin==0xffffffffu only happens to the first record
				insert_offset(idx->index[save_tid], save_bin, save_off, last_off);
			save_off = last_off;
			save_bin = last_bin = intv.bin;
			save_tid = intv.tid;
			if (save_tid < 0) break;
		}
		if (bgzf_tell(fp) <= last_off) {
			fprintf(stderr, "[ti_index_core] bug in BGZF: %llx < %llx\n",
					(unsigned long long)bgzf_tell(fp), (unsigned long long)last_off);
			exit(1);
		}
		last_off = bgzf_tell(fp);
		last_coor = intv.beg;
	}
	if (save_tid >= 0) insert_offset(idx->index[save_tid], save_bin, save_off, bgzf_tell(fp));
	merge_chunks(idx);
	fill_missing(idx);
	if (offset0 != (uint64_t)-1 && idx->n && idx->index2[0].offset) {
		int i, beg = offset0>>32, end = offset0&0xffffffffu;
		for (i = beg; i <= end; ++i) idx->index2[0].offset[i] = 0;
	}

	free(str->s); free(str);
	return idx;
}

void ti_index_destroy(ti_index_t *idx)
{
	khint_t k;
	int i;
	if (idx == 0) return;
	// destroy the name hash table
	for (k = kh_begin(idx->tname); k != kh_end(idx->tname); ++k) {
		if (kh_exist(idx->tname, k))
			free((char*)kh_key(idx->tname, k));
	}
	kh_destroy(s, idx->tname);
	// destroy the binning index
	for (i = 0; i < idx->n; ++i) {
		khash_t(i) *index = idx->index[i];
		ti_lidx_t *index2 = idx->index2 + i;
		for (k = kh_begin(index); k != kh_end(index); ++k) {
			if (kh_exist(index, k))
				free(kh_value(index, k).list);
		}
		kh_destroy(i, index);
		free(index2->offset);
	}
	free(idx->index);
	// destroy the linear index
	free(idx->index2);
	free(idx);
}

/******************
 * index file I/O *
 ******************/

void ti_index_save(const ti_index_t *idx, BGZF *fp)
{
	int32_t i, size, ti_is_be;
	khint_t k;
	ti_is_be = bam_is_big_endian();
	bgzf_write(fp, "TBI\1", 4);
	if (ti_is_be) {
		uint32_t x = idx->n;
		bgzf_write(fp, bam_swap_endian_4p(&x), 4);
	} else bgzf_write(fp, &idx->n, 4);
	assert(sizeof(ti_conf_t) == 24);
	if (ti_is_be) { // write ti_conf_t;
		uint32_t x[6];
		memcpy(x, &idx->conf, 24);
		for (i = 0; i < 6; ++i) bgzf_write(fp, bam_swap_endian_4p(&x[i]), 4);
	} else bgzf_write(fp, &idx->conf, sizeof(ti_conf_t));
	{ // write target names
		char **name;
		int32_t l = 0;
		name = calloc(kh_size(idx->tname), sizeof(void*));
		for (k = kh_begin(idx->tname); k != kh_end(idx->tname); ++k)
			if (kh_exist(idx->tname, k))
				name[kh_value(idx->tname, k)] = (char*)kh_key(idx->tname, k);
		for (i = 0; i < kh_size(idx->tname); ++i)
			l += strlen(name[i]) + 1;
		if (ti_is_be) bgzf_write(fp, bam_swap_endian_4p(&l), 4);
		else bgzf_write(fp, &l, 4);
		for (i = 0; i < kh_size(idx->tname); ++i)
			bgzf_write(fp, name[i], strlen(name[i]) + 1);
		free(name);
	}
	for (i = 0; i < idx->n; ++i) {
		khash_t(i) *index = idx->index[i];
		ti_lidx_t *index2 = idx->index2 + i;
		// write binning index
		size = kh_size(index);
		if (ti_is_be) { // big endian
			uint32_t x = size;
			bgzf_write(fp, bam_swap_endian_4p(&x), 4);
		} else bgzf_write(fp, &size, 4);
		for (k = kh_begin(index); k != kh_end(index); ++k) {
			if (kh_exist(index, k)) {
				ti_binlist_t *p = &kh_value(index, k);
				if (ti_is_be) { // big endian
					uint32_t x;
					x = kh_key(index, k); bgzf_write(fp, bam_swap_endian_4p(&x), 4);
					x = p->n; bgzf_write(fp, bam_swap_endian_4p(&x), 4);
					for (x = 0; (int)x < p->n; ++x) {
						bam_swap_endian_8p(&p->list[x].u);
						bam_swap_endian_8p(&p->list[x].v);
					}
					bgzf_write(fp, p->list, 16 * p->n);
					for (x = 0; (int)x < p->n; ++x) {
						bam_swap_endian_8p(&p->list[x].u);
						bam_swap_endian_8p(&p->list[x].v);
					}
				} else {
					bgzf_write(fp, &kh_key(index, k), 4);
					bgzf_write(fp, &p->n, 4);
					bgzf_write(fp, p->list, 16 * p->n);
				}
			}
		}
		// write linear index (index2)
		if (ti_is_be) {
			int x = index2->n;
			bgzf_write(fp, bam_swap_endian_4p(&x), 4);
		} else bgzf_write(fp, &index2->n, 4);
		if (ti_is_be) { // big endian
			int x;
			for (x = 0; (int)x < index2->n; ++x)
				bam_swap_endian_8p(&index2->offset[x]);
			bgzf_write(fp, index2->offset, 8 * index2->n);
			for (x = 0; (int)x < index2->n; ++x)
				bam_swap_endian_8p(&index2->offset[x]);
		} else bgzf_write(fp, index2->offset, 8 * index2->n);
	}
}

static ti_index_t *ti_index_load_core(BGZF *fp)
{
	int i, ti_is_be;
	char magic[4];
	ti_index_t *idx;
	ti_is_be = bam_is_big_endian();
	if (fp == 0) {
		fprintf(stderr, "[ti_index_load_core] fail to load index.\n");
		return 0;
	}
	bgzf_read(fp, magic, 4);
	if (strncmp(magic, "TBI\1", 4)) {
		fprintf(stderr, "[ti_index_load] wrong magic number.\n");
		return 0;
	}
	idx = (ti_index_t*)calloc(1, sizeof(ti_index_t));	
	bgzf_read(fp, &idx->n, 4);
	if (ti_is_be) bam_swap_endian_4p(&idx->n);
	idx->tname = kh_init(s);
	idx->index = (khash_t(i)**)calloc(idx->n, sizeof(void*));
	idx->index2 = (ti_lidx_t*)calloc(idx->n, sizeof(ti_lidx_t));
	// read idx->conf
	bgzf_read(fp, &idx->conf, sizeof(ti_conf_t));
	if (ti_is_be) {
		bam_swap_endian_4p(&idx->conf.preset);
		bam_swap_endian_4p(&idx->conf.sc);
		bam_swap_endian_4p(&idx->conf.bc);
		bam_swap_endian_4p(&idx->conf.ec);
		bam_swap_endian_4p(&idx->conf.meta_char);
		bam_swap_endian_4p(&idx->conf.line_skip);
	}
	{ // read target names
		int j, ret;
		kstring_t *str;
		int32_t l;
		uint8_t *buf;
		bgzf_read(fp, &l, 4);
		if (ti_is_be) bam_swap_endian_4p(&l);
		buf = calloc(l, 1);
		bgzf_read(fp, buf, l);
		str = calloc(1, sizeof(kstring_t));
		for (i = j = 0; i < l; ++i) {
			if (buf[i] == 0) {
				khint_t k = kh_put(s, idx->tname, strdup(str->s), &ret);
				kh_value(idx->tname, k) = j++;
				str->l = 0;
			} else kputc(buf[i], str);
		}
		free(str->s); free(str); free(buf);
	}
	for (i = 0; i < idx->n; ++i) {
		khash_t(i) *index;
		ti_lidx_t *index2 = idx->index2 + i;
		uint32_t key, size;
		khint_t k;
		int j, ret;
		ti_binlist_t *p;
		index = idx->index[i] = kh_init(i);
		// load binning index
		bgzf_read(fp, &size, 4);
		if (ti_is_be) bam_swap_endian_4p(&size);
		for (j = 0; j < (int)size; ++j) {
			bgzf_read(fp, &key, 4);
			if (ti_is_be) bam_swap_endian_4p(&key);
			k = kh_put(i, index, key, &ret);
			p = &kh_value(index, k);
			bgzf_read(fp, &p->n, 4);
			if (ti_is_be) bam_swap_endian_4p(&p->n);
			p->m = p->n;
			p->list = (pair64_t*)malloc(p->m * 16);
			bgzf_read(fp, p->list, 16 * p->n);
			if (ti_is_be) {
				int x;
				for (x = 0; x < p->n; ++x) {
					bam_swap_endian_8p(&p->list[x].u);
					bam_swap_endian_8p(&p->list[x].v);
				}
			}
		}
		// load linear index
		bgzf_read(fp, &index2->n, 4);
		if (ti_is_be) bam_swap_endian_4p(&index2->n);
		index2->m = index2->n;
		index2->offset = (uint64_t*)calloc(index2->m, 8);
		bgzf_read(fp, index2->offset, index2->n * 8);
		if (ti_is_be)
			for (j = 0; j < index2->n; ++j) bam_swap_endian_8p(&index2->offset[j]);
	}
	return idx;
}

ti_index_t *ti_index_load_local(const char *fnidx)
{
	BGZF *fp;
	fp = bgzf_open(fnidx, "r");
	if (fp) {
		ti_index_t *idx = ti_index_load_core(fp);
		bgzf_close(fp);
		return idx;
	} else return 0;
}

#ifdef _USE_KNETFILE
static void download_from_remote(const char *url)
{
	const int buf_size = 1 * 1024 * 1024;
	char *fn;
	FILE *fp;
	uint8_t *buf;
	knetFile *fp_remote;
	int l;
	if (strstr(url, "ftp://") != url && strstr(url, "http://") != url) return;
	l = strlen(url);
	for (fn = (char*)url + l - 1; fn >= url; --fn)
		if (*fn == '/') break;
	++fn; // fn now points to the file name
	fp_remote = knet_open(url, "r");
	if (fp_remote == 0) {
		fprintf(stderr, "[download_from_remote] fail to open remote file.\n");
		return;
	}
	if ((fp = fopen(fn, "w")) == 0) {
		fprintf(stderr, "[download_from_remote] fail to create file in the working directory.\n");
		knet_close(fp_remote);
		return;
	}
	buf = (uint8_t*)calloc(buf_size, 1);
	while ((l = knet_read(fp_remote, buf, buf_size)) != 0)
		fwrite(buf, 1, l, fp);
	free(buf);
	fclose(fp);
	knet_close(fp_remote);
}
#else
static void download_from_remote(const char *url)
{
	return;
}
#endif

static char *get_local_version(const char *fn)
{
    struct stat sbuf;
	char *fnidx = (char*)calloc(strlen(fn) + 5, 1);
	strcat(strcpy(fnidx, fn), ".tbi");
	if ((strstr(fnidx, "ftp://") == fnidx || strstr(fnidx, "http://") == fnidx)) {
		char *p, *url;
		int l = strlen(fnidx);
		for (p = fnidx + l - 1; p >= fnidx; --p)
			if (*p == '/') break;
		url = fnidx; fnidx = strdup(p + 1);
		if (stat(fnidx, &sbuf) == 0) {
			free(url);
			return fnidx;
		}
		fprintf(stderr, "[%s] downloading the index file...\n", __func__);
		download_from_remote(url);
		free(url);
	}
    if (stat(fnidx, &sbuf) == 0) return fnidx;
	free(fnidx); return 0;
}

const char **ti_seqname(const ti_index_t *idx, int *n)
{
	const char **names;
	khint_t k;
	*n = idx->n;
	names = calloc(idx->n, sizeof(void*));
	for (k = kh_begin(idx->tname); k < kh_end(idx->tname); ++k)
		if (kh_exist(idx->tname, k))
			names[kh_val(idx->tname, k)] = kh_key(idx->tname, k);
	return names;
}

ti_index_t *ti_index_load(const char *fn)
{
	ti_index_t *idx;
    char *fname = get_local_version(fn);
	if (fname == 0) return 0;
	idx = ti_index_load_local(fname);
	if (idx == 0) fprintf(stderr, "[ti_index_load] fail to load the index: %s\n", fname);
    free(fname);
	return idx;
}

int ti_index_build2(const char *fn, const ti_conf_t *conf, const char *_fnidx)
{
	char *fnidx;
	BGZF *fp, *fpidx;
	ti_index_t *idx;
	if ((fp = bgzf_open(fn, "r")) == 0) {
		fprintf(stderr, "[ti_index_build2] fail to open the file: %s\n", fn);
		return -1;
	}
	idx = ti_index_core(fp, conf);
	bgzf_close(fp);
	if (_fnidx == 0) {
		fnidx = (char*)calloc(strlen(fn) + 5, 1);
		strcpy(fnidx, fn); strcat(fnidx, ".tbi");
	} else fnidx = strdup(_fnidx);
	fpidx = bgzf_open(fnidx, "w");
	if (fpidx == 0) {
		fprintf(stderr, "[ti_index_build2] fail to create the index file.\n");
		free(fnidx);
		return -1;
	}
	ti_index_save(idx, fpidx);
	ti_index_destroy(idx);
	bgzf_close(fpidx);
	free(fnidx);
	return 0;
}

int ti_index_build(const char *fn, const ti_conf_t *conf)
{
	return ti_index_build2(fn, conf, 0);
}

/********************************************
 * parse a region in the format chr:beg-end *
 ********************************************/

int ti_get_tid(const ti_index_t *idx, const char *name)
{
	khiter_t iter;
	const khash_t(s) *h = idx->tname;
	iter = kh_get(s, h, name); /* get the tid */
	if (iter == kh_end(h)) return -1;
	return kh_value(h, iter);
}

int ti_parse_region(const ti_index_t *idx, const char *str, int *tid, int *begin, int *end)
{
	char *s, *p;
	int i, l, k;
	l = strlen(str);
	p = s = (char*)malloc(l+1);
	/* squeeze out "," */
	for (i = k = 0; i != l; ++i)
		if (str[i] != ',' && !isspace(str[i])) s[k++] = str[i];
	s[k] = 0;
	for (i = 0; i != k; ++i) if (s[i] == ':') break;
	s[i] = 0;
	if ((*tid = ti_get_tid(idx, s)) < 0) {
		free(s);
		return -1;
	}
	if (i == k) { /* dump the whole sequence */
		*begin = 0; *end = 1<<29; free(s);
		return 0;
	}
	for (p = s + i + 1; i != k; ++i) if (s[i] == '-') break;
	*begin = atoi(p);
	if (i < k) {
		p = s + i + 1;
		*end = atoi(p);
	} else *end = 1<<29;
	if (*begin > 0) --*begin;
	free(s);
	if (*begin > *end) return -1;
	return 0;
}

/*******************************
 * retrieve a specified region *
 *******************************/

#define MAX_BIN 37450 // =(8^6-1)/7+1

static inline int reg2bins(uint32_t beg, uint32_t end, uint16_t list[MAX_BIN])
{
	int i = 0, k;
	if (beg >= end) return 0;
	if (end >= 1u<<29) end = 1u<<29;
	--end;
	list[i++] = 0;
	for (k =    1 + (beg>>26); k <=    1 + (end>>26); ++k) list[i++] = k;
	for (k =    9 + (beg>>23); k <=    9 + (end>>23); ++k) list[i++] = k;
	for (k =   73 + (beg>>20); k <=   73 + (end>>20); ++k) list[i++] = k;
	for (k =  585 + (beg>>17); k <=  585 + (end>>17); ++k) list[i++] = k;
	for (k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) list[i++] = k;
	return i;
}

ti_iter_t ti_iter_first()
{
	ti_iter_t iter;
	iter = calloc(1, sizeof(struct __ti_iter_t));
	iter->from_first = 1;
	return iter;
}

ti_iter_t ti_iter_query(const ti_index_t *idx, int tid, int beg, int end)
{
	uint16_t *bins;
	int i, n_bins, n_off;
	pair64_t *off;
	khint_t k;
	khash_t(i) *index;
	uint64_t min_off;
	ti_iter_t iter = 0;

	if (beg < 0) beg = 0;
	if (end < beg) return 0;
	// initialize the iterator
	iter = calloc(1, sizeof(struct __ti_iter_t));
	iter->idx = idx; iter->tid = tid; iter->beg = beg; iter->end = end; iter->i = -1;
	// random access
	bins = (uint16_t*)calloc(MAX_BIN, 2);
	n_bins = reg2bins(beg, end, bins);
	index = idx->index[tid];
	if (idx->index2[tid].n > 0) {
		min_off = (beg>>TAD_LIDX_SHIFT >= idx->index2[tid].n)? idx->index2[tid].offset[idx->index2[tid].n-1]
			: idx->index2[tid].offset[beg>>TAD_LIDX_SHIFT];
		if (min_off == 0) { // improvement for index files built by tabix prior to 0.1.4
			int n = beg>>TAD_LIDX_SHIFT;
			if (n > idx->index2[tid].n) n = idx->index2[tid].n;
			for (i = n - 1; i >= 0; --i)
				if (idx->index2[tid].offset[i] != 0) break;
			if (i >= 0) min_off = idx->index2[tid].offset[i];
		}
	} else min_off = 0; // tabix 0.1.2 may produce such index files
	for (i = n_off = 0; i < n_bins; ++i) {
		if ((k = kh_get(i, index, bins[i])) != kh_end(index))
			n_off += kh_value(index, k).n;
	}
	if (n_off == 0) {
		free(bins); return iter;
	}
	off = (pair64_t*)calloc(n_off, 16);
	for (i = n_off = 0; i < n_bins; ++i) {
		if ((k = kh_get(i, index, bins[i])) != kh_end(index)) {
			int j;
			ti_binlist_t *p = &kh_value(index, k);
			for (j = 0; j < p->n; ++j)
				if (p->list[j].v > min_off) off[n_off++] = p->list[j];
		}
	}
	if (n_off == 0) {
		free(bins); free(off); return iter;
	}
	free(bins);
	{
		int l;
		ks_introsort(offt, n_off, off);
		// resolve completely contained adjacent blocks
		for (i = 1, l = 0; i < n_off; ++i)
			if (off[l].v < off[i].v)
				off[++l] = off[i];
		n_off = l + 1;
		// resolve overlaps between adjacent blocks; this may happen due to the merge in indexing
		for (i = 1; i < n_off; ++i)
			if (off[i-1].v >= off[i].u) off[i-1].v = off[i].u;
		{ // merge adjacent blocks
			for (i = 1, l = 0; i < n_off; ++i) {
				if (off[l].v>>16 == off[i].u>>16) off[l].v = off[i].v;
				else off[++l] = off[i];
			}
			n_off = l + 1;
		}
	}
	iter->n_off = n_off; iter->off = off;
	return iter;
}

const char *ti_iter_read(BGZF *fp, ti_iter_t iter, int *len)
{
	if (iter->finished) return 0;
	if (iter->from_first) {
		int ret;
		if ((ret = ti_readline(fp, &iter->str)) < 0) {
			iter->finished = 1;
			return 0;
		} else {
			if (len) *len = iter->str.l;
			return iter->str.s;
		}
	}
	if (iter->n_off == 0) return 0;
	while (1) {
		int ret;
		if (iter->curr_off == 0 || iter->curr_off >= iter->off[iter->i].v) { // then jump to the next chunk
			if (iter->i == iter->n_off - 1) break; // no more chunks
			if (iter->i >= 0) assert(iter->curr_off == iter->off[iter->i].v); // otherwise bug
			if (iter->i < 0 || iter->off[iter->i].v != iter->off[iter->i+1].u) { // not adjacent chunks; then seek
				bgzf_seek(fp, iter->off[iter->i+1].u, SEEK_SET);
				iter->curr_off = bgzf_tell(fp);
			}
			++iter->i;
		}
		if ((ret = ti_readline(fp, &iter->str)) >= 0) {
			ti_intv_t intv;
			iter->curr_off = bgzf_tell(fp);
			if (iter->str.s[0] == iter->idx->conf.meta_char) continue;
			get_intv((ti_index_t*)iter->idx, &iter->str, &intv);
			if (intv.tid != iter->tid || intv.beg >= iter->end) break; // no need to proceed
			else if (intv.end > iter->beg && iter->end > intv.beg) {
				if (len) *len = iter->str.l;
				return iter->str.s;
			}
		} else break; // end of file
	}
	iter->finished = 1;
	return 0;
}

void ti_iter_destroy(ti_iter_t iter)
{
	if (iter) {
		free(iter->str.s); free(iter->off);
		free(iter);
	}
}

int ti_fetch(BGZF *fp, const ti_index_t *idx, int tid, int beg, int end, void *data, ti_fetch_f func)
{
	ti_iter_t iter;
	const char *s;
	int len;
	iter = ti_iter_query(idx, tid, beg, end);
	while ((s = ti_iter_read(fp, iter, &len)) != 0)
		func(len, s, data);
	ti_iter_destroy(iter);
	return 0;
}

const ti_conf_t *ti_get_conf(ti_index_t *idx) { return idx? &idx->conf : 0; }

/*******************
 * High-level APIs *
 *******************/

tabix_t *ti_open(const char *fn, const char *fnidx)
{
	tabix_t *t;
	BGZF *fp;
	if ((fp = bgzf_open(fn, "r")) == 0) return 0;
	t = calloc(1, sizeof(tabix_t));
	t->fn = strdup(fn);
	if (fnidx) t->fnidx = strdup(fnidx);
	t->fp = fp;
	return t;
}

void ti_close(tabix_t *t)
{
	if (t) {
		bgzf_close(t->fp);
		if (t->idx) ti_index_destroy(t->idx);
		free(t->fn); free(t->fnidx);
		free(t);
	}
}

int ti_lazy_index_load(tabix_t *t)
{
	if (t->idx == 0) { // load index
		if (t->fnidx) t->idx = ti_index_load_local(t->fnidx);
		else t->idx = ti_index_load(t->fn);
		if (t->idx == 0) return -1; // fail to load index
	}
	return 0;
}

ti_iter_t ti_queryi(tabix_t *t, int tid, int beg, int end)
{
	if (tid < 0) return ti_iter_first();
	if (ti_lazy_index_load(t) != 0) return 0;
	return ti_iter_query(t->idx, tid, beg, end);	
}

ti_iter_t ti_querys(tabix_t *t, const char *reg)
{
	int tid, beg, end;
	if (reg == 0) return ti_iter_first();
	if (ti_lazy_index_load(t) != 0) return 0;
	if (ti_parse_region(t->idx, reg, &tid, &beg, &end) < 0) return 0;
	return ti_iter_query(t->idx, tid, beg, end);
}

ti_iter_t ti_query(tabix_t *t, const char *name, int beg, int end)
{
	int tid;
	if (name == 0) return ti_iter_first();
	// then need to load the index
	if (ti_lazy_index_load(t) != 0) return 0;
	if ((tid = ti_get_tid(t->idx, name)) < 0) return 0;
	return ti_iter_query(t->idx, tid, beg, end);
}

const char *ti_read(tabix_t *t, ti_iter_t iter, int *len)
{
	return ti_iter_read(t->fp, iter, len);
}
