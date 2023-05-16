#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <htslib/sam.h>

// read bam file and output cell gene UMI for each
void *write_cell_gene_UMI(
    char *bam_file,
    char *output_file)
{


    FILE *output = fopen(output_file, "w");

    if (output == NULL)
    {
        printf("ERROR: Cannot open output file %s\n", output_file);
        exit(1);
    }

    // open bam file
    samFile *bam_reader = hts_open(bam_file, "r");

    if (bam_reader == NULL)
    {
        printf("ERROR: Cannot open bam file %s\n", bam_file);
        exit(1);
    }

    // read bam header
    bam_hdr_t *bam_header = sam_hdr_read(bam_reader);

    // initialize bam record to stroe each read
    bam1_t *bam_record = bam_init1();

    // initialize cell barcode tree

    size_t total_read_counts = 0;

    printf("Start processing reads\n");

    // read bam file
    while (sam_read1(bam_reader, bam_header, bam_record) >= 0)
    {
        total_read_counts++;

        // printf("total_read_counts: %lu\n", total_read_counts);

        if (total_read_counts % 1000000 == 0)
        {
            // time stamp
            time_t t = time(NULL);

            // print time stamp
            struct tm *tm = localtime(&t);
            char s[64];
            strftime(s, sizeof(s), "%c", tm);
            printf("%s: ", s);

            printf("Processed %lu reads \n", total_read_counts);
        }

        uint8_t *cb = bam_aux_get(bam_record, "CB");

        // as some read may not be aligned to any gene, we need to check if TX tag exists
        if (cb != NULL)
        {
            char *cell_barcode = bam_aux2Z(cb);
            // printf("%s\n", cell_barcode);

            uint8_t *xf = bam_aux_get(bam_record, "xf");
            int UMI_quality = bam_aux2i(xf);

            if (UMI_quality == 25)
            {
                uint8_t *gx = bam_aux_get(bam_record, "GX");
                uint8_t *ub = bam_aux_get(bam_record, "UB");
                char *gene_ID = bam_aux2Z(gx);
                char *UMI = bam_aux2Z(ub);

                // write output to file
                fprintf(output, "%s,%s,%s\n", cell_barcode, gene_ID, UMI);

            }
        }
    }

    printf("Processed all %lu reads\n", total_read_counts);

    bam_destroy1(bam_record);
    bam_hdr_destroy(bam_header);
    hts_close(bam_reader);

    fclose(output);

}

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        printf("Usage: ./convert_bam <bam_file> <output_file>\n");
        exit(1);
    }

    char *bam_file = argv[1];
    char *output_file = argv[2];

    write_cell_gene_UMI(bam_file, output_file);

    return 0;
}
