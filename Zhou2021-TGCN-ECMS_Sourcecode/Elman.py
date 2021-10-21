
# Elman network
import numpy as np
#import pytorch
from tr import excel_to_matrix
from tr import unnormal
# import csv
import pandas as pd

def sigmoid(x):
    return np.tanh(x)
    #return np.maximum(0, x)


def dsigmoid(x):
    return 1.0 - x ** 2
    #return 1

class Elman:
    def __init__(self, *args):
        '''初始化给定尺寸的感知器'''
        self.shape = args
        n = len(args)

        # Build layers
        self.layers = []

        # 输入层（+1单位偏移量 + 第一个隐藏层的大小）
        # np.ones（size）-返回一个1的数组
        self.layers.append(np.ones(self.shape[0] + 1 + self.shape[1]))

        # 隐藏层+输出层
        for i in range(1, n):
            self.layers.append(np.ones(self.shape[i]))

        # 构造权重矩阵
        self.weights = []
        for i in range(n - 1):
            # 将返回矩阵
            self.weights.append(np.zeros((self.layers[i].size, self.layers[i + 1].size)))

        # dw会改变导数的最后一个变化（冲动）
        # [0，] * number-将返回一个数组，该数组的元素数乘以一个数值，然后乘以第一个因子
        self.dw = [0, ] * len(self.weights)

        # Lose weight
        self.reset()

    def reset(self):
        ''' Lose weight '''
        for i in range(len(self.weights)):
            Z = np.random.random((self.layers[i].size, self.layers[i + 1].size))
            self.weights[i][...] = (2 * Z - 1) * 0.25

    def propagate_forward(self, data):
        ''' 将数据从输入层传播到输出层 '''

        # 设置输入数据层
        self.layers[0][:self.shape[0]] = data
        # и первый скрытый слой
        # 和第一个隐藏层
        self.layers[0][self.shape[0]:-1] = self.layers[1]

        # 使用Sigmoid作为激活函数从第0层传播到第n-1层
        for i in range(1, len(self.shape)):
            # Propagate activity
             self.layers[i][...] = sigmoid(np.dot(self.layers[i - 1], self.weights[i - 1]))
        # 返回输出
        return self.layers[-1]

    def propagate_backward(self, target, lrate=0.05, momentum=0.1):
        ''' Backpropagation error related to target use of speed. '''
        deltas = []

        # Compute error on output layer
        error = target - self.layers[-1]
        delta = error * dsigmoid(self.layers[-1])
        deltas.append(delta)

        # Calculation error on hidden layers
        for i in range(len(self.shape) - 2, 0, -1):
            delta = np.dot(deltas[0], self.weights[i].T) * dsigmoid(self.layers[i])
            deltas.insert(0, delta)

        # Updating the balance
        for i in range(len(self.weights)):
            layer = np.atleast_2d(self.layers[i])
            delta = np.atleast_2d(deltas[i])
            dw = np.dot(layer.T, delta)
            self.weights[i] += lrate * dw + momentum * self.dw[i]
            self.dw[i] = dw

        # Return error
        return (error ** 2).sum()


# -----------------------------------------------------------------------------
if __name__ == '__main__':
    import matplotlib
    import matplotlib.pyplot as plt
    # -------------------------------------------------------------------------
    # 创建网络
    network = Elman(29, 30, 30, 30, 1)
    samples = np.zeros(2700, dtype=[('input', float, 29), ('output', float, 1)])
    samples1 = np.zeros(300, dtype=[('input', float, 29), ('output', float, 1)])



    datafile = u'data\\io_train.xlsx'
    datafile1 = u'data\\io_test.xlsx'
    sample_W = excel_to_matrix(datafile)
    sample1_W = excel_to_matrix(datafile1)
    train = []
    train_w = []
    test = []
    test_w = []
    output_pred = []
    output_groundtruth = []
    for i in range(2700):
        train.append(sample_W[i][0:29])
        train_w.append(sample_W[i][-1])
    for i in range(2700):
        samples[i] = train[i], train_w[i]
        # print(samples[i])
    for i in range(300):
        test.append(sample1_W[i][0:29])
        test_w.append(sample1_W[i][-1])
    for i in range(300):
        samples1[i] = train[i], train_w[i]
        # print(samples[i])
    for i in range(100000): #初始5w次
        n = i % samples.size
        network.propagate_forward(samples['input'][n])
        network.propagate_backward(samples['output'][n])
    for i in range(samples1.size):
        o = network.propagate_forward(samples['input'][i])
        output_pred.append(o.item())
        output_groundtruth.append(samples['output'][i])
        # print('Sample %d: %s -> %s' % (i, samples['input'][i], samples['output'][i]))
        # print('Network output: %s' % o.astype(float))
        # print()

##################################################################
    # output_mat = np.zeros((300,30))
    # #之前的归一化之前数据
    # # for i in range(300):
    # #     output_list.append([output_pred[i],output_groundtruth[i]])
    #
    # # print(output_list)
    # for i in range(300):
    #     n = i % samples.size
    #     output_mat[i][0] = np.mat(output_pred[i])[:np.newaxis]
    # print(output_mat)
######################################################################

    output_catmat = np.column_stack((samples1['input'], output_pred))
    output_catgroundth = np.column_stack((samples1['input'], samples1['output']))
    # print(output_catmat)
    outpput_mat = unnormal(output_catmat)
    outpput_groundth =unnormal(output_catgroundth)
    # outpput_mat2 = unnormal(output_list2)
    #     # np.mat([output_pred,output_groundtruth])
    # for i in range(300):
    #     print(outpput_mat[i][30])

    df = outpput_mat[:,-1]
    df_data = pd.DataFrame(df)
    outputfile = "io_predict.xlsx"
    df_data.to_excel(outputfile)
    for i in range(300):
        print(outpput_mat[i][-1],outpput_groundth[i][-1])