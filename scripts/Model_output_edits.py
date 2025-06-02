
import os
directory = os.getcwd() + "/Data/Model_output/regression/regression"


print("Files in '%s':" % directory)

for i in os.listdir(directory):
    
  file_path = os.path.join(directory, i)

  directory = os.getcwd() + "/Data/Model_output/regression/regression"
  print("Files in '%s':" % directory)
  for i in os.listdir(directory):
    
    file_path = os.path.join(directory, i)

    # The followng code will check the date of the file and move it to a new directory 
    # Check if the current item in the directory is a file
    if os.path.isfile(file_path):
    # Check if the file name contains "2024-06-"
    if "2024-06-" in i:
      # Extract the two-digit number after "2024-06-"
      xx = i.split("2024-06-")[1][:2]
      print(xx)
      # Check if the extracted number is a digit and equal to "08"
      if xx.isdigit() and xx == "08":
      # Create a new directory called "new_folder" if it doesn't exist
      new_directory = os.path.join(directory, "new_folder")
      if not os.path.exists(new_directory):
        os.makedirs(new_directory)
      # Move the file to the new directory
      new_file_path = os.path.join(new_directory, i)
      os.rename(file_path, new_file_path)
















        