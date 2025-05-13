import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')

def find_steady_points(array, output):

    # 计算相邻元素的差异（增量）
    data = np.cumsum(array)
    diff = np.abs(np.diff(data))
    
    threshold = np.median(diff)
    # 找到增量小于阈值的点
    steady_indices = np.where(diff < threshold)[0] + 1  # 索引从1开始
    
    # 如果找到了平稳点，取第一个平稳点的索引
    if len(steady_indices) > 0:
        first_steady_index = steady_indices[0]
    else:
        first_steady_index = None

    # 绘图
    #plt.plot(data, label='Cumulative Data')
    #plt.scatter(steady_indices, data[steady_indices], color='red', label='Steady Points')
    
    #if first_steady_index is not None:
    #    plt.axvline(x=first_steady_index, color='green', linestyle='--', label='First Steady Point')
    
    #plt.xlabel('Index')
    #plt.ylabel('Value')
    #plt.title('Steady Points in Cumulative Data')
    #plt.legend()
    #plt.savefig(f"{output}", dpi=300, bbox_inches='tight')
    #plt.close()
    
    return steady_indices[0]


class Find_Cluster_thr:
    def __init__(self, expr_df):
        self.expr_df = expr_df

    def find_by_distance(self, th):
        if th == 'auto':
            inx = find_steady_points(self.expr_df["Number of clusters"], "")
            best_th = self.expr_df.iloc[inx]['Distance threshold']
            return best_th
        else:
            return float(th)

    def find_by_cluster(self, th):
        if th == 'auto':
            raise('cluster 暂时不支持auto, 建议设定成和亚群数量一致的cluster')
        else:
            th = int(th)
            if not self.expr_df[self.expr_df['Number of clusters'] == th].empty:
                return self.expr_df[self.expr_df['Number of clusters'] == th]['Distance threshold'].to_list()[0]
            
            elif not self.expr_df[self.expr_df['Number of clusters'] == (th + 1)].empty:
                return self.expr_df[self.expr_df['Number of clusters'] == (th + 1)]['Distance threshold'].to_list()[0]
            
            elif not self.expr_df[self.expr_df['Number of clusters'] == (th -1)].empty:
                return self.expr_df[self.expr_df['Number of clusters'] == (th - 1)]['Distance threshold'].to_list()[0]
            
            elif not self.expr_df[self.expr_df['Number of clusters'] == (th + 2)].empty:
                return self.expr_df[self.expr_df['Number of clusters'] == (th + 2)]['Distance threshold'].to_list()[0]
            
            elif not self.expr_df[self.expr_df['Number of clusters'] == (th -2)].empty:
                return self.expr_df[self.expr_df['Number of clusters'] == (th - 2)]['Distance threshold'].to_list()[0]

            else:
                print(f'error cluster th: {th}')
                th = self.find_by_distance('auto')
                return th