from decimal import *
from math import log10

# ставим глобаный пресижен для работы с floating point numbers и открываем файл с нуклеотидными последовательностями
getcontext().prec = 60000

PATH_TO_DATA_FILE = 'data/hg38_chr1_and_chr2.fa'  # путь к файлу с двумя хромосомами

# получаем нуклеотиды второй хромосомы 10000000-10100000
letters = ''.join(open(PATH_TO_DATA_FILE, 'r').read().split('\n')).split('>chr2')[1][10000000:10100000]

# * A C G T
# A .......
# C .......
# G .......
# T .......
matrix_helper_dict = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

# матрицы переходов для CGi и non-CGi, полученные через 1 хромосому
CpG_matrix = [
    [.188, .273, .426, .112],
    [.158, .366, .283, .191],
    [.162, .351, .365, .121],
    [.087, .366, .356, .087]
]
non_CpG_matrix = [
    [.323, .174, .248, .253],
    [.347, .263, .047, .341],
    [.284, .212, .263, .239],
    [.214, .207, .254, .325]
]


def write_answer_down(start_index, length, s_x):
    # функция для записи верных ответов в файл
    if s_x < 0 or length < 200:
        return
    res_file = open('result.txt', 'a')
    res_file.write(
        'index: ' + str(10000000 + start_index) +
        ', S(x): ' + str(s_x) +
        ', length: ' + str(length) + '\n'
    )


def get_matrix_chance(first_letter, second_letter, matrix) -> Decimal:
    # Возвращает шанс перехода для двух нуклеотидов в определенной матрице переходов
    first_index = matrix_helper_dict[first_letter.upper()]
    second_index = matrix_helper_dict[second_letter.upper()]
    return Decimal(matrix[first_index][second_index])


def calculate_chance(letter_string: str, is_cpg: bool, base: Decimal = 1) -> Decimal:
    # Вспомогательная функция для поиска соотношения правдоподобия
    chance = Decimal(base)
    pair_index = 0
    while pair_index + 1 < len(letter_string):
        first_letter = letter_string[pair_index:pair_index + 2][0]
        second_letter = letter_string[pair_index:pair_index + 2][1]
        chance *= get_matrix_chance(first_letter, second_letter, CpG_matrix if is_cpg else non_CpG_matrix)
        if second_letter.upper() == 'N' or first_letter.upper() == 'N':
            chance = 1
        pair_index += 1
    return chance


# первый цикл, в котором мы ищем CpGi
while True:
    index = letters.find('CG')
    current_chain_length = 100
    peak_s = 0
    end = False
    last_data = (None, None)

    # Второй цикл, в котором мы ищем их продолжительность
    while True:
        island_letters = letters[index:index + current_chain_length].upper()
        c_count = island_letters.count('C')
        g_count = island_letters.count('G')
        cg_count = island_letters.count('CG')
        expected_to_real_cpg_relation = cg_count * current_chain_length / c_count / g_count

        if not (g_count + c_count > current_chain_length * .55 and expected_to_real_cpg_relation > 0.65):
            if current_chain_length > 100:
                current_chain_length -= 8
                end = True
            else:
                break

        # Оптимизация - вместо подсчета всего шанса с нуля каждый раз, если мы можем - мы используем последний шанс и
        #  умножаем его на шансы новых нуклеотидов
        if last_data == (None, None):
            cpg_chance = calculate_chance(island_letters, is_cpg=True)
            non_cpg_chance = calculate_chance(island_letters, is_cpg=False)
        else:
            cpg_chance = calculate_chance(island_letters[-9:], base=last_data[0], is_cpg=True)
            non_cpg_chance = calculate_chance(island_letters[-9:], base=last_data[1], is_cpg=False)
        last_data = (cpg_chance, non_cpg_chance)

        # s - соотношения правдоподобия
        s = log10(cpg_chance / non_cpg_chance)
        if s < 0:
            if current_chain_length > 100:
                # Сломалось не с первой попытки, т.е. раньше были успешные попытки
                current_chain_length -= 8
                end = True
            else:
                # Сломалось с первой попытки, игнорируем
                index += current_chain_length
                break

        # для дебага, чтобы было видно что программа что-то делает
        # print(index, s, current_chain_length)

        if peak_s * .7 > s and peak_s != 0:
            # S(x) упало на 30% или больше относительно максимума (для этого островка), значит островок закончился
            write_answer_down(index, current_chain_length, s)
            index += current_chain_length
            break

        if end:
            # По какому-то из критериев островок закончился
            write_answer_down(index, current_chain_length, s)
            index += current_chain_length
            break
        else:
            # Островок еще не закончился, идем дальше
            peak_s = s if s > peak_s else peak_s
            current_chain_length += 8
            continue

    # Ищем следующий возможный островок
    index = letters.find('CG', index + 1)
    if index == -1:
        # Таких возможных островков нет, заканчиваем программу
        break

