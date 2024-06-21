""" Routines to convert files and objects from a format to another """


def convert_gbk_to_faa(gbk_input_filename, faa_output_filename):
    """
    A function to convert a .gbk file returned from the RAST annotation to
    a .faa that can be used with ModelSEEDpy
    """
    from Bio import SeqIO

    input_handle    = open(gbk_input_filename, "r")
    output_handle = open(faa_output_filename, "w")

    for seq_record in SeqIO.parse(input_handle, "genbank") :
        print("Dealing with GenBank record %s" % seq_record.id)
        for entry in seq_record.features:
            if entry.type == "CDS":

                cds_dict = dict(entry.qualifiers)

                output_handle.write(">%s %s\n%s\n" % (
                    cds_dict["db_xref"][0],
                    cds_dict["product"][0],
                    cds_dict["translation"][0],)
                )

    output_handle.close()
    input_handle.close()

    return 0


