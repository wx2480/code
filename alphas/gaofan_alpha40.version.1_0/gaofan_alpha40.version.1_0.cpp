#include "include/qlSim.h"
namespace
{
class Gavin_Alpha_4th_1 : public AlphaBase
{
    public:
    explicit Gavin_Alpha_4th_1(XMLCONFIG::Element *cfg):
        AlphaBase(cfg), 
        close(dr.getData<QL_MATRIX<QL_FLOAT>>("adj_close")),   /*get open data*/
        vwap(dr.getData<QL_MATRIX<QL_FLOAT>>("adj_vwap")),
        volume(dr.getData<QL_MATRIX<QL_FLOAT>>("volume")), 
        spread_std(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.spread_std")),  /*get volume*/
        spread_mean(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.spread_mean")), 
        turnover_mean(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.turnover_mean")), 
        ask_mean(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.ask_mean")),
        bid_mean(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.bid_mean")),
        tasm(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.total_ask_size_mean")),
        tbsm(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.total_bid_size_mean")),
        last_open(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.last_open")),
        last_close(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.last_close")),
        sum_volume_price(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.sum_volume_price")),
        ask_ave_price(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.ask_ave_price")),
        bid_ave_price(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.bid_ave_price")),
        volume_above_mid(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.volume_above_mid")),
        ret(dr.getData<QL_MATRIX<QL_FLOAT>>("adj_return")),
        cap(dr.getData<QL_MATRIX<QL_FLOAT>>("cap")),
        buy_value_exlarge_order(dr.getData<QL_MATRIX<QL_FLOAT>>("wind.MoneyFlow.buy_value_exlarge_order")),
        sell_value_exlarge_order(dr.getData<QL_MATRIX<QL_FLOAT>>("wind.MoneyFlow.sell_value_exlarge_order")),
        s_mfd_inflowvolume(dr.getData<QL_MATRIX<QL_FLOAT>>("wind.MoneyFlow.s_mfd_inflowvolume")),
        Num1(cfg->getAttributeIntDefault("para1",5)),
        Num2(cfg->getAttributeIntDefault("para2",10)),
        Num3(cfg->getAttributeIntDefault("para3",70))
    {
        para_1.resize(GLOBAL::Dates.size());
        para_2.resize(GLOBAL::Dates.size());
        para_3.resize(GLOBAL::Dates.size());
        para_4.resize(GLOBAL::Dates.size());
        para_5.resize(GLOBAL::Dates.size());
        para_6.resize(GLOBAL::Dates.size());
        para_7.resize(GLOBAL::Dates.size());
        para_8.resize(GLOBAL::Dates.size());
        para_9.resize(GLOBAL::Dates.size());
        para_10.resize(GLOBAL::Dates.size());
        para_11.resize(GLOBAL::Dates.size());
        para_12.resize(GLOBAL::Dates.size());
        para_13.resize(GLOBAL::Dates.size());
        para_14.resize(GLOBAL::Dates.size());
        para_15.resize(GLOBAL::Dates.size());
        para_16.resize(GLOBAL::Dates.size());
        para_17.resize(GLOBAL::Dates.size());
        EMA_12.resize(GLOBAL::Dates.size());
        EMA_26.resize(GLOBAL::Dates.size());
        DIF.resize(GLOBAL::Dates.size());
        MACD.resize(GLOBAL::Dates.size());
        bar.resize(GLOBAL::Dates.size());
        MAVOL_5.resize(GLOBAL::Dates.size());
        MAVOL_10.resize(GLOBAL::Dates.size());
        MAVOL_5_10_delta.resize(GLOBAL::Dates.size());
        MB.resize(GLOBAL::Dates.size());
        MD.resize(GLOBAL::Dates.size());
        UP.resize(GLOBAL::Dates.size());
        DN.resize(GLOBAL::Dates.size());
        base_alpha.resize(GLOBAL::Dates.size());
        for(int di = 0; di < GLOBAL::Dates.size(); ++di)
        {
            para_1[di].resize(GLOBAL::Instruments.size(), nan);
            para_2[di].resize(GLOBAL::Instruments.size(), nan);
            para_3[di].resize(GLOBAL::Instruments.size(), nan);
            para_4[di].resize(GLOBAL::Instruments.size(), nan);
            para_5[di].resize(GLOBAL::Instruments.size(), nan);
            para_6[di].resize(GLOBAL::Instruments.size(), nan);
            para_7[di].resize(GLOBAL::Instruments.size(), nan);
            para_8[di].resize(GLOBAL::Instruments.size(), nan);
            para_9[di].resize(GLOBAL::Instruments.size(), nan);
            para_10[di].resize(GLOBAL::Instruments.size(), nan);
            para_11[di].resize(GLOBAL::Instruments.size(), nan);
            para_12[di].resize(GLOBAL::Instruments.size(), nan);
            para_13[di].resize(GLOBAL::Instruments.size(), nan);
            para_14[di].resize(GLOBAL::Instruments.size(), nan);
            para_15[di].resize(GLOBAL::Instruments.size(), nan);
            para_16[di].resize(GLOBAL::Instruments.size(), nan);
            para_17[di].resize(GLOBAL::Instruments.size(), nan);
            EMA_12[di].resize(GLOBAL::Instruments.size(), nan);
            EMA_26[di].resize(GLOBAL::Instruments.size(), nan);
            DIF[di].resize(GLOBAL::Instruments.size(), nan);
            MACD[di].resize(GLOBAL::Instruments.size(), nan);
            bar[di].resize(GLOBAL::Instruments.size(), nan);
            MAVOL_5[di].resize(GLOBAL::Instruments.size(), nan);
            MAVOL_10[di].resize(GLOBAL::Instruments.size(), nan);
            MAVOL_5_10_delta[di].resize(GLOBAL::Instruments.size(), nan);
            MB[di].resize(GLOBAL::Instruments.size(), nan);
            MD[di].resize(GLOBAL::Instruments.size(), nan);
            UP[di].resize(GLOBAL::Instruments.size(), nan);
            DN[di].resize(GLOBAL::Instruments.size(), nan);
            base_alpha[di].resize(GLOBAL::Instruments.size(), nan);
        }
    }




    

    void generate(int di) override
    {
        QL_Oputils::qilin_ar_init(para_1, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(para_2, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(para_3, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(para_4, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(para_5, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(para_6, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(para_7, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(para_8, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(para_9, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(para_10, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(para_11, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(para_12, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(para_13, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(para_14, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(para_15, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(para_16, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(para_17, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(EMA_12, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(EMA_26, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(DIF, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(MACD, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(bar, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(MAVOL_5, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(MAVOL_10, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(MAVOL_5_10_delta, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(MB, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(MD, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(UP, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(DN, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(base_alpha, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);

        calc_EMA(di);
        calc_MACD(di);
        rank_DIF_MACD(di);
        calc_MAVOL(di);
        daily_calc(di);
        cal_BOLL(di);
        cal_para(di);
        base_calc(di);

        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
        {
            if((valid[ di ][ ii ]))
            {
                std::vector<float> past_alpha(Num3,nan);
                for(int i = 0; i < Num3; i++)
                {
                    past_alpha[i] = base_alpha[di-delay-i][ii];
                }
                alpha[ii] = ts_rank(past_alpha) ;
                //alpha[ii] = base_alpha[di-delay][ii] ;

            }
        }
        return;
    }

    void base_calc(int di)
    {
        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
        {
            base_alpha[di-delay][ii] =   para_1[di-delay][ii] * para_2[di-delay][ii]  * para_13[di-delay][ii];
        }
        weighted_group_rank(di-delay, 90 ,base_alpha[di-delay]);
        return;
    }

    void daily_calc(int di)
    {
        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
        {
            std::vector<float> past_MAVOL_5_10(5,nan);
            for(int i = 0; i < 5; ++ i)
            {
                past_MAVOL_5_10[i] = MAVOL_5_10_delta[di-delay-5+i][ii]  ;

            }
            if(ask_mean(di-delay,45,ii)!=0)
            {
                para_1[di-delay][ii] = -(ask_mean(di-delay,48,ii) - ask_mean(di-delay,45,ii))/ask_mean(di-delay,45,ii) ;
            }
            para_5[di-delay][ii] = tasm(di-delay,48,ii) - tbsm(di-delay,48,ii) - (tasm(di-delay,45,ii) - tbsm(di-delay,45,ii));
            if(close[di-delay][ii]!=0)
            {
                para_2[di-delay][ii] = (vwap[di-delay][ii] - close[di-delay][ii])/ close[di-delay][ii] ;
            }
            para_3[di-delay][ii] = ts_rank(past_MAVOL_5_10) * bar[di-delay][ii];
            if(last_close(di-delay,48,ii)!=0)
            {
                para_6[di-delay][ii] = (last_open(di-delay,43,ii) - last_close(di-delay,48,ii))/last_close(di-delay,48,ii);
            }
            if(last_close(di-delay,1,ii)!=0)
            {
                para_7[di-delay][ii] = (last_open(di-delay,5,ii) - last_close(di-delay,1,ii))/last_close(di-delay,1,ii);
            }
            if(sum_volume_price(di-delay,45,ii) != 0)
            {
                para_8[di-delay][ii] = -(sum_volume_price(di-delay,48,ii) - sum_volume_price(di-delay,45,ii))/sum_volume_price(di-delay,45,ii);
            }
            if(bid_ave_price(di-delay,45,ii) != 0)
            {
                para_9[di-delay][ii] = -(bid_ave_price(di-delay,48,ii) - bid_ave_price(di-delay,45,ii))/bid_ave_price(di-delay,45,ii);
            }
            std::vector<float> past_volume_above_mid(6,nan);
            for(int i = 0; i < 6; i++)
            {
                past_volume_above_mid[i] = volume_above_mid(di-delay,42+i,ii);
            }
            para_10[di-delay][ii] = QL_Oputils::slop(past_volume_above_mid.begin(),past_volume_above_mid.end());
            if(tbsm(di-delay,48,ii)!=0)
            {
                para_11[di-delay][ii] = (tbsm(di-delay,43,ii) - tbsm(di-delay,48,ii))/tbsm(di-delay,48,ii) ;
            }
            if(ask_mean(di-delay,48,ii)!=0)
            {
                para_12[di-delay][ii] = (ask_mean(di-delay,43,ii) - ask_mean(di-delay,48,ii))/ask_mean(di-delay,48,ii) ;
            }
            para_13[di-delay][ii] = s_mfd_inflowvolume[di-delay][ii] ;
    
            

        }
        QL_Oputils::rank(para_1[di-delay]);
        QL_Oputils::rank(para_2[di-delay]);
        QL_Oputils::rank(para_3[di-delay]);
        QL_Oputils::rank(para_5[di-delay]);
        QL_Oputils::rank(para_6[di-delay]);
        QL_Oputils::rank(para_7[di-delay]);
        QL_Oputils::rank(para_8[di-delay]);
        QL_Oputils::rank(para_9[di-delay]);
        QL_Oputils::rank(para_10[di-delay]);
        QL_Oputils::rank(para_11[di-delay]);
        QL_Oputils::rank(para_12[di-delay]);
        QL_Oputils::rank(para_13[di-delay]);
        QL_Oputils::rank(para_14[di-delay]);
        QL_Oputils::rank(para_15[di-delay]);
        QL_Oputils::rank(para_16[di-delay]);
        QL_Oputils::rank(para_17[di-delay]);
        return;
    }


     void calc_EMA(int di)
    {
        int DI1 = QLCalendar::getDi(20090106);
        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
        {
            EMA_12[DI1-delay-1][ii] = 0;
            EMA_26[DI1-delay-1][ii] = 0;
        }
        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
        {
            if(isnormal(close[di-delay][ii])==1)
            {
                EMA_12[di-delay][ii] = EMA_12[di-delay-1][ii] * 11/13.0 + close[di-delay][ii] * 2/13.0;
                EMA_26[di-delay][ii] = EMA_26[di-delay-1][ii] * 25/27.0 + close[di-delay][ii] * 2/27.0;
                DIF[di-delay][ii] = EMA_12[di-delay][ii] - EMA_26[di-delay][ii];
            }
            else
            {
                EMA_12[di-delay][ii] = EMA_12[di-delay-1][ii];
                EMA_26[di-delay][ii] = EMA_26[di-delay-1][ii];
                DIF[di-delay][ii] = DIF[di-delay-1][ii];
            }
        }
        return;
    }

    void calc_MACD(int di)
    {
        int DI = QLCalendar::getDi(20090106);
        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
        {
            MACD[DI-delay-1][ii] = 0;
        }
        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
        {
            MACD[di-delay][ii] = MACD[di-delay-1][ii] * 8/10.0 + DIF[di-delay][ii] * 2/10.0;
        }
        return;
    }

    void rank_DIF_MACD(int di)
    {
        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
        {
            bar[di-delay][ii] = DIF[di-delay][ii] - MACD[di-delay][ii];        
        }
        QL_Oputils::rank(bar[di-delay]);
        return;
    }

    void calc_MAVOL(int di)
    {
        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
        {
            std::vector<float> past_5_volume(Num1,nan);
            std::vector<float> past_10_volume(Num2,nan);
            for(int i = 0; i < Num1; ++i)
            {
                past_5_volume[i] = volume[di-delay-i][ii];
            }
            for(int i = 0; i < Num2; ++i)
            {
                past_10_volume[i] = volume[di-delay-i][ii];
            }
            MAVOL_5[di-delay][ii] = QL_Oputils::mean(past_5_volume);
            MAVOL_10[di-delay][ii] = QL_Oputils::mean(past_10_volume);
        }
        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
        {
            MAVOL_5_10_delta[di-delay][ii] = (MAVOL_5[di-delay][ii] - MAVOL_10[di-delay][ii]) ;
        }
        QL_Oputils::rank(MAVOL_5_10_delta[di-delay]);
        return;
    }

    float ts_rank(vector<float> &vect)
    {
        QL_Oputils::rank(vect);
        return vect[0];
    }

    void cal_BOLL(int di)
    {
        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
        {
            std::vector<float> past_close_19(25,nan);
            std::vector<float> past_close_20(26,nan);
            for(int i = 0; i < 25; ++i)
            {
                past_close_19[i] = close[di-delay-25+i][ii];
            }
            for(int i = 0; i < 26; ++i)
            {
                past_close_20[i] = close[di-delay-26+i][ii];
            }
            MB[di-delay][ii] = QL_Oputils::mean(past_close_19);
            MD[di-delay][ii] = QL_Oputils::std(past_close_20);
            UP[di-delay][ii] = MB[di-delay][ii] + 2 * MD[di-delay][ii];
            DN[di-delay][ii] = MB[di-delay][ii] - 2 * MD[di-delay][ii];
        }
        return;
    }


    void cal_para(int di)
    {
        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
        {
            std::vector<float> past_MB_10(10,nan);
            std::vector<float> past_UP_10(10,nan);
            std::vector<float> past_DN_10(10,nan);
            for(int i = 0; i < 10; ++i)
            {
                past_MB_10[i] = MB[di-delay-10+i][ii];
                past_UP_10[i] = UP[di-delay-10+i][ii];
                past_DN_10[i] = DN[di-delay-10+i][ii];
            }
            std::vector<float> past_MAVOL_5_10(5,nan);
            for(int i = 0; i < 5; ++ i)
            {
                past_MAVOL_5_10[i] = MAVOL_5_10_delta[di-delay-5+i][ii]  ;

            }
            para_4[di-delay][ii] = ts_rank(past_MAVOL_5_10)
            * ( 1 - QL_Oputils::corr(past_UP_10,past_DN_10)) * (2 + slop_ols(past_UP_10) - slop_ols(past_DN_10));
        }
        QL_Oputils::rank(para_4[di-delay]);
        return;
    }

    float slop_ols(vector<float>&vect)
    {
        int n = vect.size();
        float mu = QL_Oputils::mean(vect);

        vector<float> vect1(n, nan);
        for(int i = 0; i < n; ++i)
        {
            vect1[i] = powf((i-(n-1)/2.0), 2);
        }
        float L_xx = QL_Oputils::sum(vect1);

        vector<float> vect2(n, nan);
        for(int i = 0; i < n; ++i)
        {
            vect2[i] = (i-(n-1)/2.0) * (vect[i]-mu);
        }
        float L_xy = QL_Oputils::sum(vect2);

        if (fabs(L_xx) > 1e-3) return L_xy / L_xx;
        else return nan;
    } 

    void weighted_group_rank(int di, int group_total, vector<float>&cross_sectional_vect)
    {
        vector<float> criterion(GLOBAL::Instruments.size(),nan);

        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ii)
        {
            criterion[ii] = cap[di-delay][ii] ;
        }
        QL_Oputils::rank(criterion);

        vector<vector<float>> grouped_data(group_total, vector<float>(GLOBAL::Instruments.size(),nan));

        vector<vector<float>> weight_standards(group_total, vector<float>(GLOBAL::Instruments.size(),nan));

        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ii)
        {
            if(!GLOBALFUNC::iserr(criterion[ii]))
            {
                int group_id = floor(criterion[ii] * group_total);

                if(group_id == group_total) group_id = group_total-1;

                grouped_data[group_id][ii] = cross_sectional_vect[ii];

                weight_standards[group_id][ii] = ret[di-delay][ii];
            }
        }

        std::vector<float> group_weight(group_total, nan);

        for(int group_id = 0; group_id < group_total; ++ group_id)
        {
            QL_Oputils::rank(grouped_data[group_id]);

            group_weight[group_id] = QL_Oputils::mean(weight_standards[group_id]);

        }
        QL_Oputils::rank(group_weight);

        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ii)
        {
            if(!GLOBALFUNC::iserr(criterion[ii]))
            {
                int group_id = floor(criterion[ii] * group_total);

                if(group_id == group_total) group_id = group_total-1;

                cross_sectional_vect[ii] = grouped_data[group_id][ii] * group_weight[group_id];
            }
        }
        return;

    }




    void checkPointSave(boost::archive::binary_oarchive &ar) 
    {
        ar & *this;
    }
    void checkPointLoad(boost::archive::binary_iarchive &ar)
    {
        ar & *this;
    }
    std::string version() const
    {
        return GLOBALCONST::VERSION;
    }
    private:
        friend class boost::serialization::access;
        template<typename Archive>
        void serialize(Archive & ar, const unsigned int/*file_version*/)
        {
            ar & para_1;
            ar & para_2;
            ar & para_3;
            ar & para_4;
            ar & para_5;
            ar & para_6;
            ar & para_7;
            ar & para_8;
            ar & para_9;
            ar & para_10;
            ar & para_11;
            ar & para_12;
            ar & para_13;
            ar & para_14;
            ar & para_15;
            ar & para_16;
            ar & para_17;
            ar & EMA_12;
            ar & EMA_26;
            ar & DIF;
            ar & MACD;
            ar & bar;
            ar & MAVOL_5;
            ar & MAVOL_10;
            ar & MAVOL_5_10_delta;
            ar & MB;
            ar & MD;
            ar & UP;
            ar & DN;
            ar & base_alpha;
        }

        vector<vector<float>> para_1;
        vector<vector<float>> para_2;
        vector<vector<float>> para_3;
        vector<vector<float>> para_4;
        vector<vector<float>> para_5;
        vector<vector<float>> para_6;
        vector<vector<float>> para_7;
        vector<vector<float>> para_8;
        vector<vector<float>> para_9;
        vector<vector<float>> para_10;
        vector<vector<float>> para_11;
        vector<vector<float>> para_12;
        vector<vector<float>> para_13;
        vector<vector<float>> para_14;
        vector<vector<float>> para_15;
        vector<vector<float>> para_16;
        vector<vector<float>> para_17;
        vector<vector<float>> EMA_12;
        vector<vector<float>> EMA_26;
        vector<vector<float>> DIF;
        vector<vector<float>> MACD;
        vector<vector<float>> bar;
        vector<vector<float>> MAVOL_5;
        vector<vector<float>> MAVOL_10;
        vector<vector<float>> MAVOL_5_10_delta;
        vector<vector<float>> MB;
        vector<vector<float>> MD;
        vector<vector<float>> UP;
        vector<vector<float>> DN;
        vector<vector<float>> base_alpha;

        const QL_MATRIX<QL_FLOAT> &close;   /*declare the data*/
        const QL_MATRIX<QL_FLOAT> &vwap;
        const QL_MATRIX<QL_FLOAT> &volume;   /*declare the data*/
        const QL_CUBE<QL_FLOAT> &spread_std;
        const QL_CUBE<QL_FLOAT> &spread_mean;
        const QL_CUBE<QL_FLOAT> &turnover_mean;
        const QL_CUBE<QL_FLOAT> &ask_mean;
        const QL_CUBE<QL_FLOAT> &bid_mean;
        const QL_CUBE<QL_FLOAT> &tasm;
        const QL_CUBE<QL_FLOAT> &tbsm;
        const QL_CUBE<QL_FLOAT> &last_open;
        const QL_CUBE<QL_FLOAT> &last_close;
        const QL_CUBE<QL_FLOAT> &sum_volume_price;
        const QL_CUBE<QL_FLOAT> &ask_ave_price;
        const QL_CUBE<QL_FLOAT> &bid_ave_price;
        const QL_CUBE<QL_FLOAT> &volume_above_mid;
        const QL_MATRIX<QL_FLOAT> &ret;
        const QL_MATRIX<QL_FLOAT> &cap;
        const QL_MATRIX<QL_FLOAT> &buy_value_exlarge_order;
        const QL_MATRIX<QL_FLOAT> &sell_value_exlarge_order;
        const QL_MATRIX<QL_FLOAT> &s_mfd_inflowvolume;
        int Num1;
        int Num2;
        int Num3;
};
}
extern "C"
{
    AlphaBase * createStrategy(XMLCONFIG::Element *cfg)
    {
        AlphaBase * str = new Gavin_Alpha_4th_1(cfg);
        return str;
    }
}

