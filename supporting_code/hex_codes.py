#generating n number of distinct hex codes for our barcodes

import random

number_of_colors = 3000 #change value as desired
color = ["#" + ''.join([random.choice('0123456789ABCDEF') for x in range(6)])
             for c in range(number_of_barcodes)]

#output: list of hex codes
