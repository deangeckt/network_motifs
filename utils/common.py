def sort_dict_freq(d: dict) -> dict:
    return dict(sorted(d.items(), key=lambda item: item[1], reverse=True))


def get_decimal_from_bin_vec(vec: list[int]) -> int:
    decimal = 0
    for i, bit in enumerate(vec):
        decimal += bit * (2 ** i)
    return decimal


def get_bin_vec_from_decimal(decimal: int, pad_to: int) -> list[int]:
    bin_digits = [int(d) for d in str(bin(decimal))[2:]]
    pad_amount = pad_to - len(bin_digits)
    padding = pad_amount * [0]
    bin_digits = padding + bin_digits
    return bin_digits
