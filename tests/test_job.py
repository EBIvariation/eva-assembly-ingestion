from eva_assembly_ingestion.job import pretty_print


def test_pretty_print():
    pretty_print(['Header 1', 'Long Header 2'],
                 [['row1 cell 1', 'row1 cell 2'], ['row2 cell 1', 'Super long row2 cell 2']])
