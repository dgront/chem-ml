class Encoder:
    def __init__(self, *allowed_labels, **kwargs):
        self.__default_idx = kwargs.get("default_idx", -1)

        self.__labels = {label: i for i, label in enumerate(allowed_labels)}
    
    def __str__(self):
        out = ""
        for label, idx in self.__labels.items():
            out += f"{label} -> {idx}\n"
        return out

    def encode_binary(self, *tokens):
        out = [0 for i in self.__labels]
        for t in tokens:
            out[self.__labels.get(t, self.__default_idx)] = 1
        return out

    def encode_values(self, *tokens):
        out = [0 for i in self.__labels]
        token_values = {label:tokens.count(label) for label in set(tokens)}
        for t, val in token_values.items(): out[self.__labels.get(t, self.__default_idx)] = val
        return out


if __name__ == "__main__":
    e = Encoder("One", "Two", "Three", "Five", default_idx=0)
    print(e)
    tokens = ["Two", "Three", "Four", "Two"]
    print(e.encode_binary(*tokens))
    print(e.encode_values(*tokens))
