import os, shutil, sys

# 为了精简输出大小，应尽可能的对文件拆分
src_dir = "./src/"
target_dir = "./scripts/pre-algo/"


# 读文件内容
def read_file(file):
    with open(file, encoding="UTF-8") as f:
        read_all = f.read()
        f.close()

    return read_all

# 写内容到文件
def rewrite_file(file, data):
    with open(file, "w", encoding="UTF-8") as f:
        f.write(data)
        f.close()

if os.path.exists(target_dir):
	shutil.rmtree(target_dir)
shutil.copytree(src_dir, target_dir)

magic = '////=======//// '

for root, ds, fs in os.walk(target_dir):
	for f in fs:
		fullname = os.path.join(root, f)
		if not fullname.endswith('.hpp'):
			continue
		content = read_file(fullname).splitlines()
		for i in range(0, len(content)):
			line = content[i]
			if line.startswith('#include <'):
				content[i] = magic + line
		content = '\n'.join(content)
		rewrite_file(fullname, content)
		sys.stderr.write(fullname + '\n')
