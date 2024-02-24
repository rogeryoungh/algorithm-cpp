import os, shutil, sys, datetime

# 为了精简输出大小，应尽可能的对文件拆分
src_path = sys.argv[1]

date = f"// GENERATE DATE: {datetime.datetime.now()}"

# 读文件内容
def read_file(file) -> str:
  with open(file, encoding="UTF-8") as f:
    read_all = f.read()
    f.close()

  return read_all

# 写内容到文件
def rewrite_file(file, data):
  with open(file, "w", encoding="UTF-8") as f:
    f.write(data)
    f.close()

def trim_include(content):
  lines = content.splitlines()
  n = len(lines)
  i = 0
  mp = set()
  while i < len(lines):
    lines[i] = lines[i].rstrip()
    if lines[i].startswith("#include <") and lines[i] in mp:
      del lines[i]
    else:
      mp.add(lines[i])
      i += 1
  return "\n".join(lines)

magic = '////=======//// '

ret = os.popen(f'g++ -E -CC -nostdinc++ -nostdinc -P -I./src -I./scripts/system-header/ {src_path}')
ret = date + '\n' + ret.read().replace(magic, '')
ret = trim_include(ret)
sys.stderr.write("success: %.3f KB\n" % (len(ret) / 1024))
print(ret)
