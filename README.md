# Games101 课程作业详解
本仓库为[GAMES101: 现代计算机图形学入门](https://sites.cs.ucsb.edu/~lingqi/teaching/games101.html)课程作业的详解，包含作业原始代码，以及每项作业的完整代码实现和思路详解，欢迎参考。
# 使用方法
## 获取作业原始代码
本仓库对GAMES101的作业原始代码按照作业序号进行了整理打包，命名为`Assignment0`-`Assignment9`，每个作业项目结构如下：
```
Assignment1
    |
    |-- Code               // 包含作业源代码
    |-- Assignment1.pdf    // 作业说明PDF
```
获取作业原始代码的命令如下：
```Bash
git clone git@github.com:sqduan/Games101.git
git checkout assignment_empty
```
## 编译
进入对应作业的Code目录下，执行编译命令，例如编译Assignment1，则执行
```
cd Assignment1/Code
mkdir build && cd build
cmake ..
make
```
## 获取作业完整实现及思路详解
> [!WARNING]  
> 在你未独自完成作业前，请尽量不要查看代码完整实现及详解

使用`git checkout`命令切换至对应的作业详解：

```
git checkout assignment1    # 切换到作业1详解
```
