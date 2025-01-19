def get_chance(numbers):
    differences = [abs(numbers[i + 1] - numbers[i]) + 1 
                   for i in range(len(numbers) - 1)]
    max_difference = max(differences)

    ratios = [difference / max_difference for difference in differences]
    chance = sum(ratios) / len(ratios)

    return chance