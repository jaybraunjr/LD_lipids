import csv
import pandas as pd
import numpy as np
import fileinput


df = pd.read_csv('xbox.txt',skiprows=10)
n = df[::100]
n.to_csv('new_xbox.csv', index=None,)

