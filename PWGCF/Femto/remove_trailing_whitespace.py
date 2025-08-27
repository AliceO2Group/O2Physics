import os

# 设置你的项目根目录
project_root = "./Tasks"  # 如果脚本放在项目根目录，"./" 即可

# 定义需要处理的文件类型，可按需扩展
file_extensions = (".h", ".hpp", ".cpp", ".c", ".cc", ".py")  

def remove_trailing_whitespace(file_path):
    with open(file_path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    
    # 删除每行末尾空格，但保留换行符
    lines = [line.rstrip() + "\n" for line in lines]

    with open(file_path, "w", encoding="utf-8") as f:
        f.writelines(lines)

def process_directory(root_dir):
    for dirpath, _, filenames in os.walk(root_dir):
        for file in filenames:
            if file.endswith(file_extensions):
                file_path = os.path.join(dirpath, file)
                remove_trailing_whitespace(file_path)
                print(f"Processed: {file_path}")

if __name__ == "__main__":
    process_directory(project_root)
    print("All done! Trailing whitespaces removed.")
