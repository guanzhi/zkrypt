
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
H[0] = 0000000000000000000000000000000000000000000000000000000000000001
H[1] = 30644e72e131a029048b6e193fd841045cea24f6fd736bec231204708f703636
H[2] = 30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000000
H[3] = 0000000000000000b3c4d79d41a91758cb49c3517c4604a520cff123608fc9cb

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
