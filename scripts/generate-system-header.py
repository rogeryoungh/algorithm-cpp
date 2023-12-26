import os, shutil, sys

src_dir = './src/'
target_dir = "./scripts/system-header/"


# 读文件内容
def read_file(file):
  with open(file, encoding="UTF-8") as f:
    read_all = f.read()
    f.close()

  return read_all

# 写内容到文件
def rewrite_file(file, data): 
    dir,fi = os.path.split(file)
    if not os.path.exists(dir):
        os.makedirs(dir)
    with open(file, "w", encoding="UTF-8") as f:
        f.write(data)
        f.close()



magic = '////=======//// '

def main():
  if os.path.exists(target_dir):
    shutil.rmtree(target_dir)
  os.mkdir(target_dir)


  header = []
  for root, ds, fs in os.walk(src_dir):
    for f in fs:
      fullname = os.path.join(root, f)
      if not fullname.endswith('.hpp'):
        continue
      content = read_file(fullname).splitlines()
      for i in range(0, len(content)):
        line = content[i]
        if line.startswith('#include <') and line.endswith('>'):
          header.append(line[len('#include <'):len(line) - 1])

  header = list(set(header))

  for h in header:
    content = magic + f'#include <{h}>\n'
    rewrite_file(target_dir + h, content)

  sys.stderr.write(str(header) + '\n')

main()
