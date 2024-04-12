def sort_dict_freq(d: dict) -> dict:
    return dict(sorted(d.items(), key=lambda item: item[1], reverse=True))

