class AlgGenerator:
    def __init__(self, a, c, m, x):
        self.a, self.c, self.m, self.x = a, c, m, x
        
    def get_sequence(self, digits_count, elements_count):
        sequence = []

        divider = pow(10, digits_count)

        while elements_count:
            self.x = (self.a * self.x + self.c) % self.m
            number = self.x

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