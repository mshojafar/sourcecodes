import numpy as np
import pandas as pd
import xlrd
from sklearn import preprocessing
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
mm = MinMaxScaler()
def excel_to_matrix(path):
    table = xlrd.open_workbook(path).sheets()[0]  # 获取第一个sheet表
    row = table.nrows  # 行数
    col = table.ncols  # 列数
    datamatrix = np.zeros((row, col))  # 生成一个nrows行ncols列，且元素均为0的初始矩阵
    for x in range(col):
        cols = np.matrix(table.col_values(x))  # 把list转换为矩阵进行矩阵操作
        datamatrix[:, x] = cols  # 按列把数据存进矩阵中
    # 数据归一化
    datamatrix1 = mm.fit_transform(datamatrix)
    return datamatrix1
def unnormal(data):
    mm_data = mm.inverse_transform(data)
    return mm_data

samples = np.zeros(300, dtype=[('input', float, 29), ('output', float, 1)])


