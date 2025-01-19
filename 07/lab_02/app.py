from PyQt5.QtGui import QPalette, QColor
from PyQt6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QPushButton, QGridLayout, QTextEdit,
    QDoubleSpinBox
)

from PyQt6.QtCore import Qt
import numpy as np

import markov

MAX_SIZE = 10

class MatrixApp(QWidget):
    num_states = 1
    matrix_widgets = []

    def __init__(self):
        super().__init__()
        
        self.initUI()

    def initUI(self):
        self.setWindowTitle("Лабораторная работа №2")
        
        # Основной вертикальный layout
        self.layout = QVBoxLayout()
        self.setStyleSheet('background-color: #F8E4F9;')

        # Виджет для ввода числа состояний
        self.state_input_label = QLabel(f"Введите количество состояний (от 1 до {MAX_SIZE}):")
        self.state_input = QDoubleSpinBox()
        self.state_input.setRange(1, MAX_SIZE)
        self.state_input.setDecimals(0)
        self.state_input.setValue(1)
        
        self.layout.addWidget(self.state_input_label)
        self.layout.addWidget(self.state_input)

        # Кнопка для создания матрицы
        self.create_matrix_button = QPushButton("Создать матрицу")
        self.create_matrix_button.clicked.connect(self.create_matrix)
        self.layout.addWidget(self.create_matrix_button)

        # Layout для матрицы
        self.matrix_layout = QGridLayout()
        self.matrix_layout .setColumnStretch(0, 1)
        for i in range(1, MAX_SIZE + 1):
            self.matrix_layout .setColumnStretch(i, 5)


        self.layout.addLayout(self.matrix_layout)
        

        # Кнопки для расчета и графика
        self.calculate_button = QPushButton("Рассчитать")
        self.calculate_button.clicked.connect(self.calculate)
        self.calculate_button.setEnabled(False)
        self.layout.addWidget(self.calculate_button)

        # Поле для вывода результатов
        self.res_layout = QGridLayout()
        self.res_layout .setColumnStretch(0, 1)
        for i in range(1, MAX_SIZE + 1):
            self.res_layout .setColumnStretch(i, 5)

        self.layout.addLayout(self.res_layout)

        self.result_output = QTextEdit()
        self.result_output.setReadOnly(True)
        self.layout.addWidget(self.result_output)

        # Устанавливаем основной layout
        self.setLayout(self.layout)

    def create_matrix(self):
        # Очищаем предыдущую матрицу
        for i in reversed(range(self.matrix_layout.count())): 
            widget = self.matrix_layout.itemAt(i).widget()
            if widget is not None:
                widget.deleteLater()
        
        # Получаем количество состояний
        try:
            self.num_states = int(self.state_input.text())
            if self.num_states < 1 or self.num_states > MAX_SIZE:
                raise ValueError(f"Число должно быть от 1 до {MAX_SIZE}.")
        except ValueError as e:
            self.result_output.setText(f"Ошибка: {str(e)}")
            return
        
        # Добавляем номера столбцов
        for j in range(self.num_states):
            column_label = QLabel(f"{j + 1}")
            self.matrix_layout.addWidget(column_label, 0, j + 1, Qt.AlignmentFlag.AlignCenter)

        # Добавляем номера строк и создаем матрицу с QDoubleSpinBox
        self.matrix_widgets = []
        for i in range(self.num_states):
            row_label = QLabel(f"{i + 1}")
            self.matrix_layout.addWidget(row_label, i + 1, 0, Qt.AlignmentFlag.AlignRight)

            row_widgets = []
            for j in range(self.num_states):
                spin_box = QDoubleSpinBox()
                spin_box.setRange(-10000.0, 10000.0)
                spin_box.setDecimals(2)
                spin_box.setValue(1.0)
                spin_box.setSingleStep(0.1)
                self.matrix_layout.addWidget(spin_box, i + 1, j + 1)
                row_widgets.append(spin_box)
            self.matrix_widgets.append(row_widgets)

        # Включаем кнопку "Результаты", если матрица создана
        if len(self.matrix_widgets) > 0:
            self.calculate_button.setEnabled(True)

    def calculate(self):
        # Инициализируем пустую матрицу
        results = []
        for row in self.matrix_widgets:
            # Создаем список для текущей строки
            current_row = []
            for widget in row:
                value = widget.value()  # Получаем значение из QDoubleSpinBox
                current_row.append(value)  # Добавляем значение в текущую строку
            results.append(current_row)  # Добавляем текущую строку в матрицу
        
        try:
            res = markov.calculate_p(np.array(results))
            time = markov.calculate_time(res, np.array(results))

            result_text = "Результаты:\n"
            for i in range(len(res)):
                result_text += f"P{i + 1} = {res[i]:.2f}\tt{i + 1} = {time[i]:.2f}\n"
            self.result_output.setText(result_text)
        except Exception as e:
            self.result_output.setText(f"Нет решения ({e}).")
            return

if __name__ == '__main__':
    app = QApplication([])
    ex = MatrixApp()
    ex.resize(400, 300)
    ex.show()
    app.exec()
    #sys.exit(app.exec_())