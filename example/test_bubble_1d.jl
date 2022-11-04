
using BubbleMsh

x_ = [0.0 0.0 0.0;1.0 0.0 0.0]
x = BubbleMsh.initnodes([0.5,0.0,0.0],[0.1,0.0,0.0],9)
x = BubbleMsh.bubble(x,x_,0.1,0.01)