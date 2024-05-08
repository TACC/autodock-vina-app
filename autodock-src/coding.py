
# Initializing list 
test_list = [('Manjeet', 10), ('Akshat', 4), ('Akash', 2), ('Nikhil', 8)] 
 
# Initializing N 
N = 2
 
# printing original list 
print("The original list is : " + str(test_list)) 
 
# Get Top N elements from Records 
# Using sorted() + lambda 
res = sorted(test_list, key = lambda x: x[1], reverse = True)[:N] 
 
# printing result 
print("The top N records are : " + str(res)) 