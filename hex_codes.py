#generating n number of distinct hex codes
number_of_colors = 3000
color = ["#" + ''.join([random.choice('0123456789ABCDEF') for x in range(6)])
             for c in range(number_of_barcodes)]
