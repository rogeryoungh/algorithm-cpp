import os, shutil, sys, datetime

# 为了精简输出大小，应尽可能的对文件拆分
algo_dir = './scripts/pre-algo/'
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

magic = '////=======//// '

content = read_file(src_path)
content = content.replace('#include "../../src', '#include "./pre-algo')
rewrite_file('./scripts/pre.full.cpp', content)

ret = os.popen('g++ -E -CC -nostdinc++ -nostdinc -P -I./scripts/algo ./scripts/pre.full.cpp')
ret = date + '\n' + ret.read().replace(magic, '')
sys.stderr.write("success: %.3f KB\n" % (len(ret) / 1024))
print(ret)
