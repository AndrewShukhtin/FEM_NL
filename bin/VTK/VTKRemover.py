import os

dir_name = os.path.abspath(os.curdir)
test = os.listdir(dir_name)

for item in test:
    if item.endswith(".vtk"):
        os.remove(os.path.join(dir_name, item))


