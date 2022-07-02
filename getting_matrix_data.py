import random
from typing import Tuple, List

import pandas as pd

CHROMOSOME_FILEPATH = 'data/hg38_chr1_and_chr2.fa'
CPG_ISLAND_CHR1_DATA_FILEPATH = 'data/cpgIslandExt_chr1.tsv'

text = open(CHROMOSOME_FILEPATH).read().split('\n')
text = ''.join(text).split('>chr2')[0]

# получаем список пар начала-конца CpGi
cpg_start_end_data = pd.read_csv(CPG_ISLAND_CHR1_DATA_FILEPATH, sep='\t')
list_of_CpG_tuples: List[Tuple[int, int]] = []
best_diff = 999999
for i in range(cpg_start_end_data.shape[0] - 1):
    diff = cpg_start_end_data['chromStart'][i + 1] - cpg_start_end_data['chromEnd'][i]
    if diff < best_diff:
        best_diff = diff
        print(best_diff, cpg_start_end_data['chromEnd'][i], cpg_start_end_data['chromStart'][i])
    list_of_CpG_tuples.append((cpg_start_end_data['chromStart'][i], cpg_start_end_data['chromEnd'][i]))

# получаем список пар начала-конца non CpGi
list_of_non_CpG_tuples = []
last_end = 0
for i in list_of_CpG_tuples:
    # для ускорения работы; разкоментируйте эту строку, чтобы использовать возможных 10% данных
    # if random.randint(1, 10) == 5:
    list_of_non_CpG_tuples.append((last_end, i[0]))
    last_end = i[1]


def get_data_collection(list_of_tuples):
    # Возвращает информацию о нуклеотидах используя текст и список начала-конца островков, которые мы будем смотреть
    data_collection = {'C': {}, 'A': {}, 'T': {}, 'G': {}}
    for xd in list_of_tuples:
        start_text_index, end_text_index = xd[0], xd[1]
        text_form = text[start_text_index:end_text_index]
        index = 0
        while index < len(text_form) - 2:
            index += 1
            letters = text_form[index:index + 2]
            first_letter = letters[0]
            second_letter = letters[1]

            # если это N, то это может быть любой нуклеотид, т.е. игнорируем
            if second_letter.upper() == 'N':
                index += 1
                continue
            elif first_letter.upper() == 'N':
                continue
            if data_collection[first_letter.upper()].get(second_letter.upper()) is None:
                data_collection[first_letter.upper()][second_letter.upper()] = 1
            else:
                data_collection[first_letter.upper()][second_letter.upper()] += 1
    return data_collection


def get_chances(data_collection_dict):
    # Возвращает матрицу вероятностей, используя информацию о нуклеотидах
    chances = data_collection_dict.copy()
    for key in data_collection_dict.keys():
        into_a = data_collection_dict[key]['A']
        into_c = data_collection_dict[key]['C']
        into_t = data_collection_dict[key]['T']
        into_g = data_collection_dict[key]['G']
        into_sum = into_g + into_t + into_c + into_a
        chances[key]['A'] = into_a / into_sum
        chances[key]['C'] = into_c / into_sum
        chances[key]['T'] = into_t / into_sum
        chances[key]['G'] = into_g / into_sum

    return chances


CpG_data_collection = get_data_collection(list_of_CpG_tuples)
CpG_chances = get_chances(CpG_data_collection)
print(CpG_chances)
# non_CpG_data_collection = get_data_collection(list_of_non_CpG_tuples)
# non_CpG_chances = get_chances(non_CpG_data_collection)
# print(non_CpG_chances)

# Результаты:
# CpG_chances
# 'C': {'A': 0.15837867821467821, 'G': 0.28384284671174637, 'C': 0.366412418078071, 'T': 0.19136605699550444}
# 'A': {'G': 0.425988734782872, 'C': 0.2734641552096582, 'T': 0.11203978492452843, 'A': 0.18850732508294135}
# 'T': {'C': 0.36641323107062657, 'T': 0.1906764131644355, 'G': 0.3556883075344337, 'A': 0.08722204823050422}
# 'G': {'C': 0.3518708776890497, 'G': 0.364915501985898, 'A': 0.16215471433620685, 'T': 0.12105890598884542}

# non_CpG_chances
# 'C': {'G': 0.046790273388744336, 'T': 0.3416565462597573, 'C': 0.26385980393730907, 'A': 0.34769337641418935}
# 'A': {'C': 0.17438652549178385, 'G': 0.24884102804407315, 'T': 0.2532893430089264, 'A': 0.3234831034552166}
# 'T': {'C': 0.2075294133033855, 'T': 0.32487080380430844, 'G': 0.25362613893546065, 'A': 0.21397364395684543}
# 'G': {'G': 0.2636327509308697, 'C': 0.21240330085286113, 'T': 0.2394685833426606, 'A': 0.2844953648736086}

