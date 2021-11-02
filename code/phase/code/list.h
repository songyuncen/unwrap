#ifndef __LIST
#define __LIST
#define POS_RES     0x01   /* 1st bit */
#define NEG_RES     0x02   /* 2nd bit */
#define VISITED     0x04   /* 3rd bit */
#define ACTIVE      0x08   /* 4th bit */
#define BRANCH_CUT  0x10   /* 5th bit */
#define BORDER      0x20   /* 6th bit */
#define UNWRAPPED   0x40   /* 7th bit */
#define POSTPONED   0x80   /* 8th bit */
#define RESIDUE      (POS_RES | NEG_RES)
#define AVOID        (BRANCH_CUT | BORDER)
int GetNextOneToUnwrap(int *a, int *b, int *index_list,
                       int *num_index, int xsize, int ysize);
void UpdateList(float *qual_map, int x, int y, float val,
                float *phase, float *soln, unsigned char *bitflags,
                int xsize, int ysize, int *index_list,
                int *num_index, int ignore_code, int processed_code,
                int postponed_code, int max_list_size, 
                int dxdy_flag, float *min_qual);
void InsertList(float *soln, float val, float *qual_map,
                unsigned char *bitflags, int a, int b,
                int *index_list, int *num_index, int xsize,
                int ysize, int processed_code, int postponed_code,
                float *min_qual, int max_list_size);
#endif
