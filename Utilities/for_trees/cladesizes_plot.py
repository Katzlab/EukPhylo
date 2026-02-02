# Author: Godwin Ani, Jan 30, 2026
#Description: Produce clade sizes plot
#Use: Python3 cladesizes.py
#use: Make sure to edit filename, column to plot, number of values to show in plots, and name of output file.
import pandas as pd
import matplotlib.pyplot as plt
import re

# Read CSV
df = pd.read_csv("LACA.csv")

values = []
for row in df["Ex_he_1"].dropna():
    row = str(row).strip()
    numbers_in_brackets = re.findall(r"\((\d+)\)", row)
    nums = map(int, re.findall(r"\d+", row))
    values.extend(nums)    
values = pd.Series(values)
values = values[values >= 3]
freq = values.value_counts().sort_index()
x_pos = range(len(freq))

# Plot
plt.bar(x_pos, freq.values)
plt.xlabel("Clade sizes")
plt.ylabel("#Clades")
plt.title("Frequency of Heterolobosea clade sizes in LACA trees")
plt.xticks(x_pos, freq.index, rotation=45, fontsize=7)
max_y = freq.max()
plt.ylim(bottom=1)
ticks = [1] + list(range(5, max_y + 1, 5))
plt.yticks(ticks)


plt.savefig('LACA_Ex_he.png', dpi=200, bbox_inches='tight')
plt.show()