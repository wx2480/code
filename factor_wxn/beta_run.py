import core
import data
import datetime
from joblib import Parallel,delayed
import numpy as np
import os, sys
import pandas as pd
import tools
from tqdm import tqdm

sys.path.append(r'/data/stock/newSystemData/feature_base/structure/yuanyang/alphasystem/dataMgr/')
import dataMgr as dM

sys.path.append(r'/data/stock/newSystemData/feature_base/structure/yili/db_utilise/')
import trading_date as ts


# **********************************************************
# *********************  based on snapshot  ****************
# **********************************************************
# Calculation Factors and save to './rawFactor/factor_name/'
class CalculationFactor:
    """
    calculate factor based on snapshot
    """
    def __init__(self, start, end, factor_list):
        self.start = start
        self.end = end
        self.factor_list = factor_list
    
    def pipiline(self, func, stkdata):
        factor = func(stkdata)
        _factor = tools.tools.de_extre(factor)
        _factor = tools.tools.normalization(_factor)
        return(factor, _factor, func)
        # MAX_WORKERS = 12
        # res = Parallel(n_jobs=MAX_WORKERS)(delayed(load_one)(path) for path in tqdm(paths))

    def run(self):
        """
        start:int
        end:int
        factor_list:a list contains some functions
        """
        timeline = ts.get_trading_date(self.start, self.end)

        for date in tqdm(timeline):
            # ********************  load data and preprocessing  ********************
            start = datetime.datetime.now()
            x1, x3, x4 = datetime.timedelta(0),datetime.timedelta(0),datetime.timedelta(0)

            stkdata = data.data._get_snap(date)
            x1 = datetime.datetime.now() - start

            time1 = datetime.datetime.now()

            save_list0 = []    # save raw factor
            save_list1 = []    # save nor factor
            
            MAX_WORKERS = 12
            res = Parallel(n_jobs=MAX_WORKERS)(delayed(self.pipiline)(func, stkdata) for func in tqdm(self.factor_list))

            # ******************** save raw data and normalization data  ********************
            for factor, _factor, func in res:
                save_list0.append((factor, date, func.__name__))
                save_list1.append((_factor, date, func.__name__))

            time2 = datetime.datetime.now()
            x3 += time2 - time1
            time5 = datetime.datetime.now()

            tools.tools.save2raw(save_list0)
            tools.tools.save2nor(save_list1)

            end = datetime.datetime.now()
            x4 = end - time5
            
            print('snapshot data')
            print('load data:  ', x1)
            print('cal factor:  ', x3)
            print('save factor:  ', x4)
            print('all time:  ', end - start)
    '''
    def run(self):
        """
        start:int
        end:int
        factor_list:a list contains some functions
        """
        timeline = ts.get_trading_date(self.start, self.end)

        for date in tqdm(timeline):
            # ********************  load data and preprocessing  ********************
            start = datetime.datetime.now()
            time1 = datetime.datetime.now()
            x1, x2, x3, x4 = datetime.timedelta(0),datetime.timedelta(0),datetime.timedelta(0),datetime.timedelta(0)

            self.stkdata = data.data.get_snap(date)
            x1 = datetime.datetime.now() - start

            save_list0 = []    # save raw factor
            save_list1 = []    # save nor factor
            for func in self.factor_list:
                time2 = datetime.datetime.now()
                factor = func(self.stkdata)
                time3 = datetime.datetime.now()
                x2 += time3 - time2

                # ******************** process factor  ********************
                _factor = tools.tools.de_extre(factor)

                _factor = tools.tools.normalization(factor)

                # ******************** save raw data and normalization data  ********************
                save_list0.append((factor, date, func.__name__))
                save_list1.append((_factor, date, func.__name__))
                time4 = datetime.datetime.now()
                x3 += time4 - time3
            time5 = datetime.datetime.now()

            tools.tools.save2raw(save_list0)
            tools.tools.save2nor(save_list1)

            end = datetime.datetime.now()
            x4 = end - time5

            print('load data:  ', x1)
            print('cal factor:  ', x2)
            print('process factor:  ', x3)
            print('save factor:  ', x4)
            print('all time:  ', end - start)
    '''
# **********************************************************
# ******************  based on transaction  ****************
# **********************************************************
# Calculation Factors and save to './rawFactor/factor_name/'
class CalculationTransactionFactor:
    def __init__(self, start, end, factor_list):
        self.start = start
        self.end = end
        self.factor_list = factor_list
    
    def pipiline(self, func, stkdata):
        factor = func(stkdata)
        _factor = tools.tools.de_extre(factor)
        _factor = tools.tools.normalization(_factor)
        return(factor, _factor, func)
        # MAX_WORKERS = 12
        # res = Parallel(n_jobs=MAX_WORKERS)(delayed(load_one)(path) for path in tqdm(paths))

    def run(self):
        """
        start:int
        end:int
        factor_list:a list contains some functions
        """
        timeline = ts.get_trading_date(self.start, self.end)

        for date in tqdm(timeline):
            # ********************  load data and preprocessing  ********************
            start = datetime.datetime.now()
            x1, x3, x4 = datetime.timedelta(0),datetime.timedelta(0),datetime.timedelta(0)

            stkdata = data.data._get_transaction(date)
            x1 = datetime.datetime.now() - start

            time1 = datetime.datetime.now()

            save_list0 = []    # save raw factor
            save_list1 = []    # save nor factor
            
            MAX_WORKERS = 12
            res = Parallel(n_jobs=MAX_WORKERS)(delayed(self.pipiline)(func, stkdata) for func in tqdm(self.factor_list))

            # ******************** save raw data and normalization data  ********************
            for factor, _factor, func in res:
                save_list0.append((factor, date, func.__name__))
                save_list1.append((_factor, date, func.__name__))

            time2 = datetime.datetime.now()
            x3 += time2 - time1
            time5 = datetime.datetime.now()

            tools.tools.save2raw(save_list0)
            tools.tools.save2nor(save_list1)

            end = datetime.datetime.now()
            x4 = end - time5

            print('transaction data')
            print('load data:  ', x1)
            print('cal factor:  ', x3)
            print('save factor:  ', x4)
            print('all time:  ', end - start)


if __name__ == "__main__":
    a = datetime.datetime.now()
    
    start = 20190603
    end = 20200824
        
    stk_list_transaction = [
        core.liquid.bar_skew_second,
        core.liquid.bar_volatility,
        core.longshort.max_buy_order,
        core.longshort.max_sell_order,
        core.longshort.amount_order,
        core.longshort.amount_buy_order,
        core.longshort.amount_sell_order,
        core.longshort.maxvol_price,
        core.volbar.volatility,
        core.volbar.bar_skew,
        core.volbar.max_gain,
        core.volbar.max_loss,
        core.volbar.break_high,
        core.volbar.break_low,
        core.volbar.normal,
        core.volbar.cross,
        core.volbar.abs_road, # 17
    ]

    Transaction = CalculationTransactionFactor(
        start,
        end,
        stk_list_transaction,
    )
    print('transaction start')
    Transaction.run()
    
    # ********************************************************************
    '''
    start = 20190603
    end = 20200814
    
    stk_list_snap =  [
        core.longshort.daily_maxdrawdown,
        core.longshort.daily_meandrawdown,
        core.longshort.daily_maxgain,
        core.longshort.daily_meangain,
        core.longshort.break_highpoint,
        core.longshort.break_lowpoint,
    ]
    
    SnapShot = CalculationFactor(
        start,
        end,
        stk_list_snap,
        )
    print('snapshot start')
    SnapShot.run()
    '''
    # ********************************************************
    print('\nit started {}.'.format(a))
    print('\nit ended {}\n'.format(datetime.datetime.now()))
