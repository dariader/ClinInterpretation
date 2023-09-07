from BamClass import BamClass


class CnvClass(BamClass):
    def __init__(self, args, bam_file):
        super().__init__(args)
        self.bam_file = bam_file

    def select_border_reads(self):
        pass

    def create_softclip_dict(self):
        pass

    def append_ends_border_reads(self):
        pass

    def append_starts_border_reads(self):
        pass

    def write_file(self):
        pass