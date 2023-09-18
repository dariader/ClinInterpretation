from BamClass import BamClass
import pysam
import numpy as np

class CnvClass(BamClass):
    """
    For read in region we need to edit:
    1) read sequence (attach clips of appropriate length from the clip dictionary) (Make reads on border paired)
    2) check the read.flag?? is should not be 'unmapped'
    3) cigar string - parse and add S counts
    4) TLEN string - ength of the reference segment mapped by the two paired-end reads.

    Action order:
    1) find reads on the region borders, create lef and right clip dictionary [done]
    2) select n of them [skip for a while]
    3) select one left read and one right (read.flag & 0x10 != 0), make them paired (exchange rnext values as start positions)
    4) attach clips to the reads (select clip size), save the position within read where the clip was inserted [done]
    5) find position in cigar string, exchange it to S<num of clip nucleotides> [done]
    """

    def __init__(self, args, bam_file):
        super().__init__(args)
        self.left_clips = []
        self.left_reads = []
        self.left_reads = []
        self.right_clips = []
        self.right_reads = []
        self.right_reads = []
        self.between_reads = []
        self.bam_file = pysam.AlignmentFile(bam_file, "rb")
        self.out_file = pysam.AlignmentFile('./TMP/temp_out.bam', "wb", template=self.bam_file)
        self.run_modification()
        self.bam_file.close()
        self.out_file.close()

    def run_modification(self):
        self.start = self.mutation_file.position.values.tolist()[0] - 1
        self.end = self.mutation_file.position_end.values.tolist()[0] - 1
        print(self.start, self.end)
        self.select_border_reads()
        left_reads = self.append_left_border_clips()
        right_reads = self.append_right_border_clips()
        for i in left_reads:
            self.out_file.write(i)
        for i in right_reads:
            self.out_file.write(i)
        for i in self.between_reads:
            if np.random.choice([0, 1], size=1, p=[0.5, 0.5]):
                self.out_file.write(i)


    def edit_cigar(self, cigar, pos, orient):
        """
         sample input - [(4, 21), (0, 69)]
        :param cigar:
        :param pos: index in the read
        :return:
        """
        new_cigar = []
        read_len = sum([i[1] for i in cigar])
        if orient == 'left':
            # [(4, 21),(2, 43), (0, 69)] --> [(4, 21),(2, 8)(4, 69+35)]
            clip_cigar = [4, 0]
            clip_len = read_len - pos
            for i in cigar:
                if clip_cigar[1] == clip_len:
                    break
                if pos == 0:
                    clip_cigar[1] += i[1]
                elif i[1] < pos: # (4, 21)
                    pos -= i[1] # 35
                    new_cigar.append(i)
                else:
                    clip_cigar[1] += i[1] - pos # 43-35
                    i = (i[0], pos) # (2, 35)
                    pos = 0
                    if i[1] > 0:
                        new_cigar.append(i)
            new_cigar.append(tuple(clip_cigar))

        else:
            # [(4, 21),(2, 43), (0, 69)] --> [(4, 56),(2, 8)(0, 69)]
            clip_cigar = [4, 0]
            print(pos)
            for i in cigar:
                if clip_cigar[1] == pos:
                    new_cigar.append(i)
                else:
                    if i[1] < pos:
                        pos -= i[1]  # 56 - 21 = 35
                        clip_cigar[1] += i[1]  # (4, 21)
                    else:
                        clip_cigar[1] += i[1]-pos  # (4, 21 + 35)
                        i = (i[0], pos)  # (2, 43-35)
                        new_cigar.append(i)
            new_cigar = [tuple(clip_cigar)] + new_cigar
        new_read_len = sum([i[1] for i in new_cigar])
        assert read_len == new_read_len, breakpoint()
        return new_cigar

    def is_right_directed(self, read):
        # Check if the 0x10 bit is set (i.e., the read is on the reverse strand)
        return read.flag & 0x10 != 0

    def select_border_reads(self):
        for read in self.bam_file:
            if self.start in read.get_reference_positions():
                pos = self.start - read.get_reference_positions()[0]
                self.left_clips.append(read.seq[:pos])
                self.left_reads.append(read)
            elif self.end in read.get_reference_positions():
                pos = self.end - read.get_reference_positions()[0]
                self.right_clips.append(read.seq[pos:])
                self.right_reads.append(read)
            elif self.start in read.get_reference_positions() and self.end in read.get_reference_positions():
                print('Long read')
                continue
            else:
                self.between_reads.append(read)

    def create_softclip_dict(self):
        pass

    def append_left_border_clips(self):
        """
        ........<---X====-----
        .......----X====>..........
        :return:
        """
        new_reads = []
        for read in self.left_reads:
            position_in_read = self.start - read.pos
            read_len = len(read.seq)
            clip_length = len(read.seq) - position_in_read  # +- 1?
            try:
                clip = [i for i in self.left_clips if len(i) >= clip_length][0][:clip_length]
            except IndexError:
                print(f'Needed clip len = {clip_length} while max clip seq = {len(max(self.left_clips))}')
                continue
            read.seq = read.seq[:position_in_read] + clip
            read.next_reference_start = position_in_read+200
            #read.seq = read.seq[:position_in_read] + 'A'*clip_length
            read.flag |= 0x100
            try:
                read.cigar = self.edit_cigar(read.cigar, position_in_read, 'left')
            except ValueError:
                breakpoint()
            assert len(read.seq) == read_len, breakpoint()
            new_reads.append(read)
        return new_reads

    def append_right_border_clips(self):
        """
         :return:
        """
        new_reads = []
        for read in self.right_reads:
            position_in_read = self.end - read.pos
            if position_in_read == 0:
                continue
            read_len = len(read.seq)
            clip_length = position_in_read # +- 1?
            try:
                clip = [i for i in self.right_clips if len(i) >= clip_length][0][:clip_length]
            except IndexError:
                breakpoint()
                print(f'Clip len = {clip_length} while read len = {read_len}')
                continue # case when there is no read of proper length
            read.seq = clip + read.seq[position_in_read:]
            read.cigar = self.edit_cigar(read.cigar, position_in_read, 'right')
            assert len(read.seq) == read_len
            new_reads.append(read)
        return new_reads

    def write_file(self):
        pass
