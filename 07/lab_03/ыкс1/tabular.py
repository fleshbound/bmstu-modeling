class TabularGenerator:
    def __init__(self, file_path):
        self.table = open(file_path, 'r')

    def get_sequence(self, digits_count, elements_count):
        sequence = []

        divider = pow(10, digits_count)

        while elements_count:
            number = self.table.read(6)

            if number == '':
                self.table.seek(0)
                number = self.table.read(6)

            number = int(number[:5])

            while number:
                if elements_count:
                    random_number = number % divider

                    if len(str(random_number)) == digits_count:
                        sequence.append(random_number)
                        elements_count -= 1

                    number //= 10
                else:
                    return sequence

        return sequence

    def __del__(self):
        self.table.close()