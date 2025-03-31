import numpy as np

def one_hot_encode(char_list, input_string):
    """
    One-hot encode a string based on a given list of characters.
    
    Parameters:
        char_list (list or str): List of unique characters defining the encoding space.
        input_string (str): The string to encode.
        
    Returns:
        np.ndarray: A 2D numpy array representing the one-hot encoded string.
    """
    # Create a mapping from character to index
    char_to_index = {char: idx for idx, char in enumerate(char_list)}
    num_chars = len(char_list)
    
    # Initialize an empty matrix
    one_hot_matrix = np.zeros((len(input_string), num_chars), dtype=int)
    
    # Fill in the one-hot matrix
    for i, char in enumerate(input_string):
        if char in char_to_index:
            one_hot_matrix[i, char_to_index[char]] = 1
        else:
            raise ValueError(f"Character '{char}' not found in the provided character list.")
    
    return one_hot_matrix

if __name__ == "__main__":
#    char_list = list("ACDEFGHIKLMNPQRSTVWY")  	# Amino acids, for instance
    char_list = list("12HUL.FIB")		# codes from PDBTM
    input_string = ".BBBBBB2BBBBBB11."

    encoded = one_hot_encode(char_list, input_string)
    print(encoded)
