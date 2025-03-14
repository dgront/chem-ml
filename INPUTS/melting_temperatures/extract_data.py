import re, sys

#
# Data extracted from:
# https://www.aatbio.com/data-sets/boiling-point-bp-and-melting-point-mp-reference-table

def extract_table_rows(html_text):
    # Pattern to match the entire <tr ...>...</tr> block with the specific attributes.
    row_pattern = r'<tr\s+style="display:none"\s+class="">(.*?)</tr>'
    # Use DOTALL so that the dot matches newline characters if needed.
    rows = re.findall(row_pattern, html_text, re.DOTALL)
    return rows

def extract_columns(row_html):
    # Pattern to extract content inside each <td>...</td> tag.
    td_pattern = r'<td>(.*?)</td>'
    # Find all <td> entries in the row.
    columns = re.findall(td_pattern, row_html, re.DOTALL)
    # Optionally, trim whitespace from each cell.
    return [col.strip() for col in columns]

def main():
    # Read the HTML file as text.
    with open(sys.argv[1], "r", encoding="utf-8") as file:
        html_text = file.read()

    # Extract each row.
    rows = extract_table_rows(html_text)
    data = []
    
    # Extract the column values for each row.
    for row in rows:
        columns = extract_columns(row)
        data.append(columns)
    
    # Print each row with its columns.
    for row in data:
        print("{}\t{}\t{}\t{}".format(*row))

if __name__ == '__main__':
    main()

