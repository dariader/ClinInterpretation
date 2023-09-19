import array
import pysam
import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("target_reads", type=str)
parser.add_argument("output_bam", type=str)
parser.add_argument("cnv_file", type=str)
arguments = parser.parse_args()
target_reads = arguments.target_reads
output_bam = arguments.output_bam
cnv_file = arguments.snp_file


# Считает покрытие в однонуклеотидной позиции
def read_cov_one_nucl(samfile, chr_name, nucl):
    cov_iter = samfile.pileup(chr_name, nucl, nucl + 1, truncate=True)
    for position in cov_iter:
        cov = position.nsegments
    return (cov)


# Сохраняет кусок в +/- 150 оснований от краев делеции
## Переписать без дублирования кода
## Сделать вменяемо с единицами
def read_del_tails(samfile, chr_name, nucl, border,
                   read_length=250):  # Добавить длину рида в инпут на старте или вычислять из бама
    const = int(read_length / 5 * 3)  # В нашем случае - 150 нуклеотидов
    if border == "left":
        reads = samfile.fetch(chr_name, nucl - const - 1, nucl - 1)
        for read in reads:
            try:
                number_list = read.get_reference_positions()
                rrs = read.get_reference_sequence()
                ind2 = number_list.index(nucl - 1)  # Потому что BAM 0-based и мы отступаем на 1 нуклетид от края
                ind1 = number_list.index(
                    nucl - const - 1)  # Потому что BAM 0-based и мы отступаем на 1 нуклетид от края, но конечный коээфициент должен быть на 1 больше
                result = rrs[ind1:ind2].upper()
                if len(result) >= 150:
                    return result
                    break
            except ValueError:
                continue
    elif border == "right":
        reads = samfile.fetch(chr_name, nucl, nucl + const)
        for read in reads:
            try:
                number_list = read.get_reference_positions()
                rrs = read.get_reference_sequence()
                ind1 = number_list.index(
                    nucl)  # Потому что BAM 0-based, но мы отступаем на 1 нуклеотид от делеции вправо
                ind2 = number_list.index(
                    nucl + const + 2)  # Потому что BAM 0-based и мы отступаем на 1 нуклеотид от делеции вправо, но конечный коээфициент должен быть на 1 больше
                result = rrs[ind1:ind2].upper()
                if len(result) >= 150:
                    return result
                    break
            except ValueError:
                continue
    else:  # Поменять на raise
        print("Border have to be only 'left' or 'right'")


# Клепаем сплит-риды
def create_split_reads(samfile_input, left_3_5, right_3_5, counter, chr_name, start, stop, read_length=250):
    # Рандомно выбираем длину концов (100-150 с каждого края)
    left_num = np.random.choice(range(100, 151), size=1)[0]
    right_num = 250 - left_num
    # Создаем сплит-рид и добавляем его свойства
    new_split_read = pysam.AlignedSegment(header=samfile_input.header)  # Надо сделать нормальный хедер
    new_split_read.query_name = f"SRRread_split_{counter}"  # ? Это норм? А для близнецового рида?
    new_split_read.query_sequence = f"{left_3_5[-left_num:]}{right_3_5[:right_num]}"
    new_split_read.flag = np.random.choice([163, 83], size=1)  # ?  Допинфа от Полины и Кати
    new_split_read.reference_name = chr_name
    if chr_name[3:] == "X":
        new_split_read.reference_id = 22
    elif chr_name[3:] == "Y":
        new_split_read.reference_id = 23
    elif type(chr_name[3:]) == int:
        new_split_read.reference_id = int(chr_name[3:]) - 1

    new_split_read.reference_start = start - 1 - left_num  # Потому что BAM 0-based
    new_split_read.mapping_quality = 20  # Все изменится, ждем данных.
    new_split_read.cigartuples = [(0, 250)]
    new_split_read.next_reference_id = new_split_read.reference_id
    new_split_read.next_reference_start = 199  # ?  Уточнить
    new_split_read.template_length = 167  # ? Уточнить
    new_split_read.query_qualities = array.array('B', np.random.choice(range(18, 36),
                                                                       size=read_length).tolist())  # ? Смущает вертикальный вывод листа
    new_split_read.tags = (("NM", 1),
                           ("RG", "L1"))  # ?   Допинфа от Полины и Кати

    # Цепляем этот же сплит-рид с другой стороны делеции
    new_split_read2 = pysam.AlignedSegment(header=samfile_input.header)  # Надо сделать нормальный хедер
    new_split_read2.query_name = new_split_read.query_name  # ? Это норм? А для близнецового рида?
    new_split_read2.query_sequence = f"{left_3_5[-left_num:]}{right_3_5[:right_num]}"
    new_split_read2.flag = new_split_read.flag  # ?  Допинфа от Полины и Кати
    new_split_read2.reference_name = chr_name
    new_split_read2.reference_id = new_split_read.reference_id
    new_split_read2.reference_start = stop - left_num  # Потому что BAM 0-based
    new_split_read2.mapping_quality = 20  # Все изменится, ждем данных.
    new_split_read2.cigar = [(0, 250)]
    new_split_read2.next_reference_id = new_split_read2.reference_id
    new_split_read2.next_reference_start = 199  # ?  Уточнить
    new_split_read2.template_length = 167  # ? Уточнить
    new_split_read2.query_qualities = new_split_read.query_qualities
    new_split_read2.tags = new_split_read.tags  # ?   Допинфа от Полины и Кати

    return new_split_read, new_split_read2


def long_insert_pair(samfile_input, chr_name, start, stop, counter,
                     min_ins_size=20, max_ins_size=200, search_radius=50, read_length=250):
    distance_to_deletion_l = np.random.choice(range(min_ins_size, max_ins_size), size=1)
    distance_to_deletion_r = np.random.choice(range(min_ins_size, max_ins_size), size=1)

    reads_l = samfile_input.fetch(chr_name, start - distance_to_deletion_l - read_length - search_radius,
                                  start - distance_to_deletion_l)

    for read_l in reads_l:
        if (read_l.reference_end < start) and (read_l.flag in [99, 163]):
            break

    reads_r = samfile_input.fetch(chr_name, stop + distance_to_deletion_r,
                                  stop + distance_to_deletion_r + search_radius + read_length)
    for read_r in reads_r:
        if (read_r.reference_start > stop) and (read_r.flag in [147, 83]):
            break

    new_read_l = pysam.AlignedSegment(header=samfile_input.header)
    new_read_r = pysam.AlignedSegment(header=samfile_input.header)
    new_read_l.query_name = f"SRR_long_insert_read_{counter}"
    new_read_l.query_sequence = read_l.query_sequence
    new_read_l.flag = read_l.flag
    new_read_l.reference_name = read_l.reference_name
    new_read_l.reference_id = read_l.reference_id
    new_read_l.reference_start = read_l.reference_start
    new_read_l.mapping_quality = read_l.mapping_quality
    new_read_l.cigarstring = read_l.cigarstring
    new_read_l.cigartuples = read_l.cigartuples
    new_read_l.next_reference_id = read_r.reference_id
    new_read_l.next_reference_start = read_r.reference_start - 1
    new_read_l.template_length = read_r.reference_end + 1 - read_l.reference_start
    new_read_l.query_qualities = read_l.query_qualities
    new_read_l.tags = read_l.tags
    new_read_r.query_name = f"SRR_long_insert_read_{counter}"
    new_read_r.query_sequence = read_r.query_sequence
    new_read_r.flag = read_r.flag
    new_read_r.reference_name = read_r.reference_name
    new_read_r.reference_id = read_r.reference_id
    new_read_r.reference_start = read_r.reference_start
    new_read_r.mapping_quality = read_r.mapping_quality
    new_read_r.cigarstring = read_r.cigarstring
    new_read_r.cigartuples = read_r.cigartuples
    new_read_r.next_reference_id = read_l.reference_id
    new_read_r.next_reference_start = read_l.reference_start - 1
    new_read_r.template_length = - (read_r.reference_end + 1 - read_l.reference_start)
    new_read_r.query_qualities = read_r.query_qualities
    new_read_r.tags = read_r.tags

    return new_read_l, new_read_r


def insertion(samfile_input, samfile_output, chr_name, start, stop):
    reads = samfile_input.fetch(chr_name, start, stop + 1)
    for read in reads:
        if read.reference_start < start:
            samfile_output.write(read)
        elif read.reference_end is None:
            if (read.reference_start + read.query_length - 1) > stop:
                samfile_output.write(read)
            else:
                if np.random.choice([0, 1], size=1) == 0:
                    samfile_output.write(read)
                else:
                    samfile_output.write(read)
                    samfile_output.write(read)
        elif read.reference_end > stop:
            samfile_output.write(read)
        else:
            if np.random.choice([0, 1], size=1) == 0:
                samfile_output.write(read)
            else:
                samfile_output.write(read)
                samfile_output.write(read)



def deletion(samfile_input, samfile_output, chr_name, start, stop, mode="hetero",
             lp_number=3):  # Добавить переменную гомо-гетеро для делеции и lp_number для количества "парочек"

    # Computing coverage, number of split reads and sequence of reference around deletion
    coverage = (int(sum(map(lambda x: read_cov_one_nucl(samfile_input, chr_name, x), [start, stop])) / 2))
    split_num = int(coverage / 5)
    left_3_5 = read_del_tails(samfile_input, chr_name, start, border="left")  # Вычисляем 3/5 длины рида до делеции
    right_3_5 = read_del_tails(samfile_input, chr_name, stop, border="right")  # Вычисляем 3/5 длины рида после делеции

    # Creating split reads
    counter = 0
    for split_read in range(split_num):
        counter += 1
        new_read, new_read2 = create_split_reads(samfile_input, left_3_5, right_3_5, counter, chr_name, start, stop)
        samfile_output.write(new_read)
        samfile_output.write(new_read2)

    # Creating pair reads with long insertion
    counter = 0
    used_list = []
    for long_pair in range(lp_number):
        # print(f"del: {used_list}")
        counter += 1
        pr1, pr2 = long_insert_pair(samfile_input, chr_name, start, stop, counter)
        samfile_output.write(pr1)
        samfile_output.write(pr2)

    # Decreasing coverage
    ## Избавиться от дублирования кода
    ## Постараться уменьшить вложенность
    reads = samfile_input.fetch(chr_name, start, stop + 1)
    for read in reads:
        if read.reference_start < start:
            samfile_output.write(read)
        elif read.reference_end is None:
            if (read.reference_start + read.query_length - 1) > stop:
                samfile_output.write(read)
            else:
                if mode == "hetero":
                    if np.random.choice([0, 1], size=1) == 0:
                        samfile_output.write(read)
                elif mode == "homo":
                    pass  # Это норм практика?
                else:  # Сделать нормальный raise
                    print("Mode of deletion may be only 'homo' or 'hetero'")
        elif read.reference_end > stop:
            samfile_output.write(read)
        else:
            if mode == "hetero":
                if np.random.choice([0, 1], size=1) == 0:
                    samfile_output.write(read)
            elif mode == "homo":
                pass  # Это норм практика?
            else:  # Сделать нормальный raise
                print("Mode of deletion may be only 'homo' or 'hetero'")


cnvs = pd.read_csv(cnv_file)

with pysam.AlignmentFile(target_reads, "rb") as samfile_input, pysam.AlignmentFile(output_bam, "wb", template=samfile_input) as samfile_output:
    chroms = list(set(samfile_input.references) & set(cnvs["chromosome"]))
    true_cnvs = cnvs.query("chromosome in @chroms")
    for cnv_index in true_cnvs.index:
        cnv_type, chr_name = true_cnvs.loc[cnv_index, "type"], true_cnvs.loc[cnv_index, "chromosome"]
        start, stop = true_cnvs.loc[cnv_index, "position_start"], true_cnvs.loc[cnv_index, "position_finish"]
        if cnv_type == "del":
            deletion(samfile_input, samfile_output, chr_name, start, stop)
        elif cnv_type == "ins":
            insertion(samfile_input, samfile_output, chr_name, start, stop)