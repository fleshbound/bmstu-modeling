from prettytable import PrettyTable

import sys

from tabular import TabularGenerator
from alg import AlgGenerator
from entropy import shannon_entropy

from PyQt5.QtWidgets import QApplication, QMainWindow, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget, QPushButton, QLineEdit, QLabel

COUNT = 10

a = 69069
c = 5
m = 2**32
start = 1

tbl_gen = TabularGenerator("./data/table.txt")
alg_gen = AlgGenerator(a, c, m, start)

# Класс основного окна
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("lab_03")
        
        self.table = QTableWidget(11, 7)  # 10 строк и 7 столбцов
        self.table.setHorizontalHeaderLabels(["Табл. 1", "Табл. 2", "Табл. 3", 
                                               "Алг. 1", "Алг. 2", "Алг. 3", 
                                               "Ручной ввод"])
        
        self.calculate_button = QPushButton("Рассчитать", self)
        self.calculate_button.clicked.connect(self.calculate_entropy)

        layout = QVBoxLayout()
        layout.addWidget(self.table)
        layout.addWidget(self.calculate_button)

        container = QWidget()
        container.setLayout(layout)
        self.setCentralWidget(container)

        self.populate_table()

    def populate_table(self):
        # Пример последовательностей (замените на свои данные)
        t_seq1 = tbl_gen.get_sequence(1, COUNT)
        t_seq2 = tbl_gen.get_sequence(2, COUNT)
        t_seq3 = tbl_gen.get_sequence(3, COUNT)
        a_seq1 = alg_gen.get_sequence(1, COUNT)
        a_seq2 = alg_gen.get_sequence(2, COUNT)
        a_seq3 = alg_gen.get_sequence(3, COUNT)

        sequences = [t_seq1, t_seq2, t_seq3, a_seq1, a_seq2, a_seq3]
        
        for i in range(len(sequences)):
            for j in range(len(sequences[i])):
                self.table.setItem(j, i, QTableWidgetItem(str(sequences[i][j])))
        
        # Рассчитываем и отображаем энтропию
        entropies = [
            f"{shannon_entropy(t_seq1):.2f}",
            f"{shannon_entropy(t_seq2):.2f}",
            f"{shannon_entropy(t_seq3):.2f}",
            f"{shannon_entropy(a_seq1):.2f}",
            f"{shannon_entropy(a_seq2):.2f}",
            f"{shannon_entropy(a_seq3):.2f}"
        ]
        
        for i in range(len(entropies)):
            self.table.setItem(10, i, QTableWidgetItem(entropies[i]))
    
    def get_user_sequence(self):
        user_sequence = []
        for row in range(self.table.rowCount() - 1):  # Проходим по всем строкам, кроме последней
            item = self.table.item(row, 6)  # Индекс 6 — это последний столбец
            if item is not None and item.text().strip():  # Проверяем, что элемент существует и не пустой
                try:
                    number = int(item.text())  # Пробуем преобразовать текст в число
                    user_sequence.append(number)
                except ValueError:
                    continue  # Игнорируем строки, которые не могут быть преобразованы в число
        return user_sequence

    def calculate_entropy(self):
        try:
            # Считываем пользовательскую последовательность из последнего столбца
            #user_sequence = list(map(int, user_sequence_text.split(',')))
            user_sequence = self.get_user_sequence()
            entropy_value = shannon_entropy(user_sequence)
            self.table.setItem(10, 6, QTableWidgetItem(f"{entropy_value:.2f}"))  # Показать в следующей строке последнего столбца
        except ValueError:
            self.table.setItem(10, 6, QTableWidgetItem("Ошибка!"))  # Показать ошибку в случае неверного ввода

# Запуск приложения
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())

'''
def main():

    tbl_gen = TabularGenerator("./data/table.txt")
    alg_gen = AlgGenerator(a, c, m, start)

    table = PrettyTable()

    #table.field_names = ["Табличный 1", "Табличный 2", "Табличный 3", "Алгоритмический 1", "Алгоритмический 2", "Алгоритмический 3"]

    t_seq1 = tbl_gen.get_sequence(1, COUNT)
    t_seq2 = tbl_gen.get_sequence(2, COUNT)
    t_seq3 = tbl_gen.get_sequence(3, COUNT)

    a_seq1 = alg_gen.get_sequence(1, COUNT)
    a_seq2 = alg_gen.get_sequence(2, COUNT)
    a_seq3 = alg_gen.get_sequence(3, COUNT)

    table.add_column("Табличный 1", t_seq1)
    table.add_column("Табличный 2", t_seq2)
    table.add_column("Табличный 3", t_seq3)

    table.add_column("Алгоритмический 1", a_seq1)
    table.add_column("Алгоритмический 2", a_seq2)
    table.add_column("Алгоритмический 3", a_seq3)
    table.add_row([
    f"{shannon_entropy(t_seq1):.2f}",
    f"{shannon_entropy(t_seq2):.2f}",
    f"{shannon_entropy(t_seq3):.2f}",
    f"{shannon_entropy(a_seq1):.2f}",
    f"{shannon_entropy(a_seq2):.2f}",
    f"{shannon_entropy(a_seq3):.2f}"
    ])
    #table.add_row([shannon_entropy(t_seq1), shannon_entropy(t_seq2), shannon_entropy(t_seq3),
    #               shannon_entropy(a_seq1), shannon_entropy(a_seq2), shannon_entropy(a_seq3)])

    print(table)

main()
'''