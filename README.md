
Zkrypt是一个开源的C语言零知识证明算法库，旨在向用户提供简洁、高效的非交互式零知识证明协议接口，用户可以通过调用接口实现完整的零知识证明协议的流程，包括公共参数设置、证明生成和验证等步骤。

本项目由北京大学关志的密码学研究组开发维护。

### 特性

- 支持多种零知识证明协议（包括Groth16, Plonk等）
- 通过算法优化提供零知识证明协议的高效实现
- 支持包括BN254、BLS381等在内的多种主流曲线，特别地支持国密SM2、SM9等算法中的推荐曲线，可与国密算法兼容
- 支持X86、ARM等硬件运行环境及Linux、Windows、Mac、Android等软件运行环境

### 零知识证明

零知识证明通常是指一种方法，其中的一个参与方（证明者）可以向另一方（验证者）证明某一论断为真（例如：拥有某一数学问题的一组解），而不泄露关于“此论断为真”以外的任何信息。注意到在此例子中，一个朴素的证明方法就是直接公开拥有的解，因此难点是在不泄露任何关于解的信息的同时实现证明。

目前已有的零知识证明协议大多具有如下形式：对于某个问题![](http://latex.codecogs.com/gif.latex?P)和值![](http://latex.codecogs.com/gif.latex?y)，证明者拥有![](http://latex.codecogs.com/gif.latex?x)使得![](http://latex.codecogs.com/gif.latex?P(x)=y)。证明者利用![](http://latex.codecogs.com/gif.latex?x)计算出一组数据（证明）并交给验证者验证，如果验证者证实这些数据确实满足协议中给出的特定关系（通常是一组等式），则验证者相信证明者确实拥有![](http://latex.codecogs.com/gif.latex?x)使得![](http://latex.codecogs.com/gif.latex?P(x)=y)。

零知识证明协议必须满足以下性质：


1. 完整性：如果![](http://latex.codecogs.com/gif.latex?P(x)=y)，则证明者使用![](http://latex.codecogs.com/gif.latex?x)生成的证明总应该被接受；
2. 可靠性：对于任何不满足![](http://latex.codecogs.com/gif.latex?P(x)=y)的![](http://latex.codecogs.com/gif.latex?x)，证明者使用![](http://latex.codecogs.com/gif.latex?x)生成的证明最多以一个小概率![](http://latex.codecogs.com/gif.latex?\rho)被接受（这使得验证者可以通过多次要求证明的方式鉴别虚假的证明者）；
3. 零知识：验证者不能从证明中获得关于![](http://latex.codecogs.com/gif.latex?x)的任何信息。

如果在一个零知识证明协议中，证明者仅需要向验证者发送证明数据，则称该协议是非交互式的（否则称为交互式的；交互式协议一般需要证明者和验证者之间传输若干随机数）。在本项目中实现的协议均为非交互式的。

### 接口功能与使用样例

此处以一个具体实例展示本算法库提供的各个接口功能与输出示例。

问题：![](http://latex.codecogs.com/gif.latex?x^2+y^2=z^2)

证明者拥有的解：![](http://latex.codecogs.com/gif.latex?(x,y,z)=(3,4,5))

使用协议：Plonk算法

#### 公共参考串生成

调用此接口可以生成零知识证明算法运行所需的必要公共参数数组。

```c
CRS[0] = 0e0a77c19a07df2f666ea36f7879462c0a78eb28f5c70b3dd35d438dc58f0d9d 1c14ef83340fbe5eccdd46def0f28c5814f1d651eb8e167ba6ba871b8b1e1b3a
CRS[1] = 17a4d6f6d16cefe72bf3edc7c37d30287a1e166c2ecb1f0715ba0e833a9e2bde 233f14020e2b1456d88258f611db0060fd1e043e1ed5ce9d065bf1b1f457d18b
CRS[2] = 16df741bb03d7bb16f1d142113831a4e2ccf1829439817cc203db9132b874414 1308727182a66078ab43f686aeb1ee4cb7f0b6468a607ae6d2b0afc861485108
CRS[3] = 06704705d990c4af35fefb49320e9d3d9f5c3235df2c1fcf3ae61e8797ba79b6 1eb00bfafd2416b9c8f3342ff00f0f29e89d295d3215b5fb9fb2de509554732c
CRS[4] = 0382c2982ed4417e5f1aa529297a93393faa918677655d94b1df0571329d18e0 20c0329323011a38975d6ce581ad8d55b4a2578035f3933b616c453ab25e9585
CRS[5] = 1eae8f756f5a8b1b6aa910784e3840e0c8832a26af578cb034327716b7c8cedd 2fc869263fb3ec9585c8250193116df69817fc4c234220c9d343109a8484c7b0
CRS[6] = 09804ebafb6af93d898c0290e4b321d20df43edb3d204815ca9e9a161ad93808 0db54595becace45889e13f1d40c6e682bceaf0b1ec76037bae22e9a4153c1db
CRS[7] = 2d61c208c72ff0cd16d553b4bca1cb4ccc27eb8d62ba1073ed2635076408d94f 230e95d188df063de791bd404b4cef67322e2cf174f0e85d195aad2df54165ad
CRS[8] = 2f0a0ced38e57517bb2092df2d74c6d5c0fbc7338aa47ded711caa39deb3a24b 287f3a9dcde3523896beb48430bad76e0a01e500c86c82aa73c7c3c9cbf59e21
CRS[9] = 0cf0cc123f00e5c888501cb41b70b3493095f818661c3744c989c2037098446d 036fff432cb0ee7c880dbff7c86b7efc81677cd5475242c4f84e361f60b859de

```

#### 证明生成

调用此接口可以让证明者将待证明的问题以电路表示，并利用其拥有的问题的一组解，根据所选用的零知识证明算法生成一个关于该问题及其解的零知识证明

```c
proof:
 A         = 015076413a1e3ccd8f8e06a3e3f3593473c7fb6c75ebee4f481465cb91f0a858 189ff80237a68c56e9a04255b09b36bfff44d0643b2551aa8a7fe7f50319bc0a
 B         = 258f280f6c3b4d89047293b190884e674350c8b5105a95b35e8500d57fb66dee 1100e2e49a1a696d0dbad87b07763f631abb256ec7613ef3112066f193b6bf93
 C         = 209cb88857f7f3691b389559bbd92e46c73f69a96ce30fed4ded56e988f2777c 10786bdf254bb11422e65a019ef465f660f00e7bb815c45c351b3bdc0b0430dc
 Z         = 056e290547f77e93184a191aea1dd4f7e8b9a6628caa3478bf93a0c138cfb72d 11b027fba02dc51a7f8d092b3d2e9c32f205538abbce60ffadbfd9a7ef6d5577
 T_lo      = 0b077787871b2a77242c4064854fcbedb2f0397ef9525f0d92afac3539dce2d5 0d99aeb6fd0f712891ec7a01f27dd514ca4a194bf93655ee6f53064747e5c83f
 T_mid     = 16ff27db58ee317a728125be457db2ae65da591663211bfed0146b7114b11c61 2bd82bfa91d559705b1fc4153bd149524106fb244fde83f444a50b730c1f1b5e
 T_hi      = 0e624a95aaa9973d0d04a41abfcd187944c65a10ef3f43f0574b05fe56c30414 08a1c1db8b10604f7980d64ad121d8b806bd407183bae466eb5769a3ccbbf519
 W_xi      = 0631d246d568f2dcda81257288be0abe5b6be4ec33b16383bdfd945f3dbf8147 0c812cd688f72b74ca5b79fcdf5bcf146f2da816134bc77be5e64de69e21b1ab
 W_xi_w    = 09d8cad01d18aae6a141abea7f4bf0d05a9f1dfca12d38f7766c53189b5f42b4 0ae73daef13721014524ac695faefdc4d523849f1391af236a87505fa0f66e17
 _a        = 21b575ce967712d9add2e169321307b03f6a99fb0beebd7dd1a63402833f5981
 _b        = 02847226776e91d204387f0d7068644b77c1e12c751907c6ffe32e54067cf549
 _c        = 013dae72fb447023df84e4bd41fbfe7efaa26fb9fc3feedd53186f58c9973c92
 _s_sigma1 = 125e521ea6c2acab9d00b1ae4bdd177ff39c5b2352e572af392ad3869ba39ed0
 _s_sigma2 = 1efb4c64bf9f43d86b583cd3ad3020718fe89c9b11aa726456dae6671eaea387
 _r        = 0f68d978fb6fb49b3dac4e86d2638e571cb1a5ced0e8ab9b9d69dd25a5868ff7
 _z_w      = 06421e1705bbbeaa3098d0afe2b85b3bc0443d1de9493b7a7b1b53b0d375490a
```

#### 证明验证

调用此接口可以让验证者验证证明者提交的证明是否正确，即证明者是否确实拥有待证明问题的一组解

```c
verify: passed
```


### 编译及使用

使用本算法库需要安装gmssl-v3或者gmssl。

#### 编译

```
gmssl@ubuntu:~/zkrypt$ mkdir build
gmssl@ubuntu:~/zkrypt$ cd build/
gmssl@ubuntu:~/zkrypt/build$ cmake ..
-- The C compiler identification is GNU 7.5.0
-- The CXX compiler identification is GNU 7.5.0
-- Check for working C compiler: /usr/bin/cc
-- Check for working C compiler: /usr/bin/cc -- works
-- Detecting C compiler ABI info
-- Detecting C compiler ABI info - done
-- Detecting C compile features
-- Detecting C compile features - done
-- Check for working CXX compiler: /usr/bin/c++
-- Check for working CXX compiler: /usr/bin/c++ -- works
-- Detecting CXX compiler ABI info
-- Detecting CXX compiler ABI info - done
-- Detecting CXX compile features
-- Detecting CXX compile features - done
-- Configuring done
-- Generating done
-- Build files have been written to: /home/gmssl/zkrypt/build
gmssl@ubuntu:~/zkrypt/build$ make
Scanning dependencies of target zkrypt
[ 12%] Building C object CMakeFiles/zkrypt.dir/src/bn254.c.o
[ 25%] Building C object CMakeFiles/zkrypt.dir/src/bn254_params.c.o
[ 37%] Building C object CMakeFiles/zkrypt.dir/src/plonk.c.o
[ 50%] Linking C shared library lib/libzkrypt.so
[ 50%] Built target zkrypt
Scanning dependencies of target plonk_test
[ 62%] Building C object CMakeFiles/plonk_test.dir/tests/plonk_test.c.o
[ 75%] Linking C executable bin/plonk_test
[ 75%] Built target plonk_test
Scanning dependencies of target bn254_test
[ 87%] Building C object CMakeFiles/bn254_test.dir/tests/plonk_test.c.o
[100%] Linking C executable bin/bn254_test
[100%] Built target bn254_test
```

#### 安装

运行make install即可安装至系统lib目录下
```
gmssl@ubuntu:~/zkrypt/build$ sudo make install
[sudo] password for gmssl: 
[ 50%] Built target zkrypt
[ 75%] Built target plonk_test
[100%] Built target bn254_test
Install the project...
-- Install configuration: ""
-- Installing: /usr/local/lib/libzkrypt.so.1.0
-- Up-to-date: /usr/local/lib/libzkrypt.so.1
-- Up-to-date: /usr/local/lib/libzkrypt.so
gmssl@ubuntu:~/zkrypt/build$ sudo ldconfig
```


#### 测试

运行make test可以查看测试样例的执行效果

```
gmssl@ubuntu:~/zkrypt/build$ make test
Running tests...
Test project /home/gmssl/zkrypt/build
    Start 1: plonk_test
1/2 Test #1: plonk_test .......................   Passed    0.20 sec
    Start 2: bn254_test
2/2 Test #2: bn254_test .......................   Passed    0.20 sec

100% tests passed, 0 tests failed out of 2

Total Test time (real) =   0.40 sec
gmssl@ubuntu:~/zkrypt/build$ 
```

