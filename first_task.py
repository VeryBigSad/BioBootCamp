from decimal import *
from math import log10

getcontext().prec = 60000

SEQUENCE_PATH = 'data/generated_sequence.fa'

nucleotides = ''.join(open(SEQUENCE_PATH, 'r').read().split('\n'))

# * A C G T
# A .......
# C .......
# G .......
# T .......
matrix_helper_dict = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
CpG_matrix = [
    [.180, .274, .426, .120],
    [.171, .368, .274, .188],
    [.161, .339, .375, .125],
    [.079, .355, .384, .182]
]
non_CpG_matrix = [
    [.300, .205, .285, .210],
    [.322, .298, .078, .302],
    [.248, .246, .298, .208],
    [.177, .239, .292, .292]
]


def get_matrix_chance(first_letter, second_letter, matrix) -> Decimal:
    # возвращает шанс из матрицы для двух нуклеотидов
    first_index = matrix_helper_dict[first_letter]
    second_index = matrix_helper_dict[second_letter]
    return Decimal(matrix[first_index][second_index])


def calculate_chance(letter_string: str, is_cpg: bool) -> Decimal:
    # вспомогательная функция при подсчете S(x)
    chance = Decimal(1)
    pair_index = 0
    while pair_index + 2 < len(letter_string):
        first_letter = letter_string[pair_index:pair_index + 2][0]
        second_letter = letter_string[pair_index:pair_index + 2][1]
        chance *= get_matrix_chance(first_letter, second_letter, CpG_matrix if is_cpg else non_CpG_matrix)
        pair_index += 1
    return chance


answer_list = []
index = nucleotides.find('CG')
while True:
    our_nucleotides = nucleotides[index:index + 300]
    g_count = our_nucleotides.count('G')
    c_count = our_nucleotides.count('C')
    if g_count + c_count > 150:
        cg_count = our_nucleotides.count('CG')
        number = cg_count * 300 / c_count / g_count
        if number > 0.6:
            # correct number
            cpg_chance = calculate_chance(our_nucleotides, is_cpg=True)
            non_cpg_chance = calculate_chance(our_nucleotides, is_cpg=False)
            s = log10(cpg_chance / non_cpg_chance)
            if s < 0:
                index += 300
                continue
            answer_list.append((index, s))
            index += 300

    index = nucleotides.find('CG', index + 1)
    if index == -1:
        break

# записываем результаты в файл
res_file = open('result.txt', 'w')
res_file.write('\n'.join([str('index: ' + str(i[0]) + ', S(x): ' + str(i[1])) for i in answer_list]))
