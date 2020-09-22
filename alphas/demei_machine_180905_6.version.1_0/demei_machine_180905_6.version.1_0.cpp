#include "include/qlSim.h"
#include <iostream>
#include <cstring>
namespace
{
class alpha03 : public AlphaBase
{
    public:
    explicit alpha03(XMLCONFIG::Element *cfg): 
        AlphaBase(cfg),
        
        
        
        member(cfg->getAttributeStringDefault("para1","123456789")),
        //Dataallbegin 
        close(dr.getData<QL_MATRIX<QL_FLOAT>>("adj_close")),
        ask_high(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.ask_high")),
        bid_high(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.bid_high")),
        last_close(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.last_close")),
        last_mean(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.last_mean")),
        mid_std(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.mid_std")),
        mid_mean(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.mid_mean")),
        spread_high(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.spread_high")),
        spread_low(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.spread_low")),
        vol_sum(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.vol_sum")),
        ask_ave_price(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.ask_ave_price")),
        bid_ave_price(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.bid_ave_price")),
        sum_volume_price(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.sum_volume_price")),
        buy_volume_small_order(dr.getData<QL_MATRIX<QL_FLOAT>>("wind.MoneyFlow.buy_volume_small_order")),
        close_net_inflow_rate_volume(dr.getData<QL_MATRIX<QL_FLOAT>>("wind.MoneyFlow.close_net_inflow_rate_volume"))
        //Dataallend

    {
        Ind.resize(5,0);
        int j =0;
        for(int i=0; i<(int)member.size();++i)
        {
            if(member.at(i) =='1')
            {
                Ind[j] = i+1;
                ++j;
            }
        }
        Alpha.resize(GLOBAL::Dates.size());
        for(int di = 0; di < GLOBAL::Dates.size(); ++di)
        {
            Alpha[di].resize(GLOBAL::Instruments.size(),nan);
            
        }
      
        
    }


    void generate(int di) override
    {
        vector<vector<float>> ALPHA;
       

        ALPHA.resize(200);

        for(int i = 0; i < 200; ++i)
        {
            ALPHA[i].resize(GLOBAL::Instruments.size(), nan);
        }
        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii) // loop each stock
        {
            if((valid[ di ][ ii ]))
            {
                ALPHA[0][ii] = 1;
            }
        }
        QL_Oputils::qilin_ar_init(Alpha, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        

        if (member.at(14-1)!='0')
        {
            for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii) // loop each stock
            {
                if((valid[ di ][ ii ]))
                {

                    ALPHA[14][ii] = 2 * close_net_inflow_rate_volume[di-delay-5][ii] - close_net_inflow_rate_volume[di-delay][ii] - close_net_inflow_rate_volume[di-delay-5*2][ii];
                }

            }
            QL_Oputils::rank(ALPHA[14]);
            for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii) // loop each stock
            {
                if((valid[ di ][ ii ]))
                {
                    if (isnan(ALPHA[14][ii])){ALPHA[14][ii] = 0.5; }
                }
            }
        }
        else
        {
            for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii) // loop each stock
            {
                if((valid[ di ][ ii ]))
                {
                    ALPHA[14][ii] = 1;
                }
            }
        }






        if (member.at(30-1) !='0')
        {
            for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii) // loop each stock
            {
                
                
                if((valid[ di ][ ii ]))
                { 
                    ALPHA[30][ii] = -(f(last_close,di-delay,ii) + f(last_close,di-delay-2,ii) - 2*f(last_close,di-delay-1,ii) );
                }
            }
            group_rank(di,20,buy_volume_small_order,ALPHA[30]);
            QL_Oputils::rank(ALPHA[30]);
            for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii) // loop each stock
            {
                if((valid[ di ][ ii ]))
                {
                    if(isnan(ALPHA[30][ii])) ALPHA[30][ii] = 0.5;
                }
            }
        }
        else
        {
            for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii) // loop each stock
            {
                if((valid[ di ][ ii ]))
                {
                    ALPHA[30][ii] = 1;
                }
            }
        }


        if (member.at(46-1) !='0')
        {
            std::vector<float> past1(12,nan);
            std::vector<float> past2(12,nan);
            std::vector<float> stat(GLOBAL::Instruments.size(),nan);
            std::vector<float> ret20day(GLOBAL::Instruments.size(),nan);
            for(int ii=0;ii<GLOBAL::Instruments.size();++ii)
            {
                
                for (int ti = 1; ti < 13; ++ti)
                {
                    past1[ti] = mid_std(di-delay,ti,ii) /(mid_mean(di-delay,ti,ii)+0.00001);
                    past2[ti] = bid_high(di-delay,ti,ii) - bid_ave_price(di-delay,ti,ii);
                }
                
                
                if(fabs(close[di-delay-15][ii])) 
                        ret20day[ii] = (close[di-delay][ii] - close[di-delay-15][ii])/close[di-delay-15][ii]; 

                stat[ii] =  -QL_Oputils::corr(past1,past2);
            }
            float Beta = QL_Oputils::beta_range(stat.begin(),stat.end(),ret20day.begin(),ret20day.end());
            for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii) // loop each stock
                {
                
                    if((valid[ di ][ ii ]))
                    {
                        ALPHA[46][ii] = stat[ii] - Beta * ret20day[ii];
                    }
                }
        }
        else
        {
            for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii) // loop each stock
            {
                if((valid[ di ][ ii ]))
                {
                    ALPHA[46][ii] = 1;
                }
            }
        }
        if (member.at(47-1) !='0')
        {
            std::vector<float> past1(12,nan);
            std::vector<float> past2(12,nan);
            std::vector<float> stat(GLOBAL::Instruments.size(),nan);
            std::vector<float> ret20day(GLOBAL::Instruments.size(),nan);
            for(int ii=0;ii<GLOBAL::Instruments.size();++ii)
            {
                
                for (int ti = 1; ti < 13; ++ti)
                {
                    past1[ti] = vol_sum(di-delay,ti,ii);
                    past2[ti] = sum_volume_price(di-delay,ti,ii)/(last_mean(di-delay,ti,ii)+0.00001);
                }
                
                
                if(fabs(close[di-delay-15][ii])) 
                        ret20day[ii] = (close[di-delay][ii] - close[di-delay-15][ii])/close[di-delay-15][ii]; 

                stat[ii] =  -QL_Oputils::corr(past1,past2);
            }
            float Beta = QL_Oputils::beta_range(stat.begin(),stat.end(),ret20day.begin(),ret20day.end());
            for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii) // loop each stock
                {
                
                    if((valid[ di ][ ii ]))
                    {
                        ALPHA[47][ii] = stat[ii] - Beta * ret20day[ii];
                    }
                }
        }
        else
        {
            for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii) // loop each stock
            {
                if((valid[ di ][ ii ]))
                {
                    ALPHA[47][ii] = 1;
                }
            }
        }
        if (member.at(51-1) !='0')
        {
            std::vector<float> past1(48,nan);
            std::vector<float> past2(48,nan);
            std::vector<float> stat(GLOBAL::Instruments.size(),nan);
            std::vector<float> ret20day(GLOBAL::Instruments.size(),nan);
            for(int ii=0;ii<GLOBAL::Instruments.size();++ii)
            {
                
                for (int ti = 1; ti < 49; ++ti)
                {
                    past1[ti] = spread_high(di-delay,ti,ii) / (spread_low(di-delay,ti,ii)+0.00001);
                    past2[ti] = ask_high(di-delay,ti,ii) - ask_ave_price(di-delay,ti,ii);
                }
                
                
                if(fabs(close[di-delay-15][ii])) 
                        ret20day[ii] = (close[di-delay][ii] - close[di-delay-15][ii])/close[di-delay-15][ii]; 

                stat[ii] =  -QL_Oputils::corr(past1,past2);
            }
            float Beta = QL_Oputils::beta_range(stat.begin(),stat.end(),ret20day.begin(),ret20day.end());
            for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii) // loop each stock
                {
                
                    if((valid[ di ][ ii ]))
                    {
                        ALPHA[51][ii] = stat[ii] - Beta * ret20day[ii];
                    }
                }
        }
        else
        {
            for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii) // loop each stock
            {
                if((valid[ di ][ ii ]))
                {
                    ALPHA[51][ii] = 1;
                }
            }
        }
        std::vector<float> power((int)member.size(),nan);
        for (int i = 0; i < (int)member.size(); ++i)
        {
            power[i] = (float)member.at(i)-48;
            // printf("%f\n", power[i]);
        }
        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii) // loop each stock
        {
            if((valid[ di ][ ii ]))
            {
                Alpha[di][ii] = (1+ALPHA[Ind[0]][ii])*(1+ALPHA[Ind[1]][ii])*(1+ALPHA[Ind[2]][ii]) ;
            }
        }
        
                
        weighted_group_rank(di,20,Alpha[di],ALPHA[Ind[3]],ALPHA[Ind[4]]);
         
        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii) // loop each stock
        {
            if((valid[ di ][ ii ]))
            {
                // std::vector<float> past_alpha(10,nan);
                // for (int j = 0; j < 10; ++j)
                // {
                //     past_alpha[j] = Alpha[di-j][ii];
                // }
                // alpha[ii] = QL_Oputils::slop(past_alpha.begin(),past_alpha.end());
                // QL_Oputils::rank(past_alpha);
                // alpha[ii] = past_alpha[0];
                alpha[ii] = Alpha[di][ii];
            }
        }
        
        return;
    }


    float f(const QL_CUBE<QL_FLOAT>& X, int di,int ii)
    {

        return(X(di,48,ii)+X(di,47,ii)+X(di,46,ii));
    }

    void group_rank(int di,int group_total, const QL_MATRIX<QL_FLOAT> &X,vector<float>&cross_sectional_vect)
    {
        vector<float> criterion(GLOBAL::Instruments.size(),nan);
        for(int ii=0;ii<GLOBAL::Instruments.size();++ii)
        {
            criterion[ii] = X[di-delay][ii];
        }
        QL_Oputils::rank(criterion);
        vector<vector<float>> grouped_data(group_total,vector<float>(GLOBAL::Instruments.size(),nan));
       
        for(int ii=0;ii<GLOBAL::Instruments.size();++ii)
        {
            if(!GLOBALFUNC::iserr(criterion[ii]))
            {
                int group_id = floor(criterion[ii]*group_total);
                if(group_id==group_total) group_id = group_total-1;
                grouped_data[group_id][ii] = cross_sectional_vect[ii];
            }
        }
        
        for(int group_id=0;group_id<group_total;++group_id)
        {
            QL_Oputils::rank(grouped_data[group_id]);
            
        }
       
        for(int ii=0;ii<GLOBAL::Instruments.size();++ii)
        {
            if(!GLOBALFUNC::iserr(criterion[ii]))
            {
                int group_id = floor(criterion[ii]*group_total);
                if(group_id==group_total) group_id = group_total-1;
                cross_sectional_vect[ii] = grouped_data[group_id][ii];
            }
        }
        return;
    }
    void weighted_group_rank(int di,int group_total,vector<float>&cross_sectional_vect,vector<float>&cri_vec,vector<float>&wei_vec)
    {
        vector<float> criterion(GLOBAL::Instruments.size(),nan);
        for(int ii=0;ii<GLOBAL::Instruments.size();++ii)
        {
            criterion[ii] = cri_vec[ii];
        }
        QL_Oputils::rank(criterion);
        vector<vector<float>> grouped_data(group_total,vector<float>(GLOBAL::Instruments.size(),nan));
        vector<vector<float>> weight_standards(group_total,vector<float>(GLOBAL::Instruments.size(),nan));
        for(int ii=0;ii<GLOBAL::Instruments.size();++ii)
        {
            if(!GLOBALFUNC::iserr(criterion[ii]))
            {
                int group_id = floor(criterion[ii]*group_total);
                if(group_id==group_total) group_id = group_total-1;
                grouped_data[group_id][ii] = cross_sectional_vect[ii];
                weight_standards[group_id][ii] = wei_vec[ii];
            }
        }
        vector<float> group_weight(group_total,nan);
        for(int group_id=0;group_id<group_total;++group_id)
        {
            QL_Oputils::rank(grouped_data[group_id]);
            group_weight[group_id] = QL_Oputils::mean(weight_standards[group_id]);
            
        }
        QL_Oputils::rank(group_weight);
        for(int ii=0;ii<GLOBAL::Instruments.size();++ii)
        {
            if(!GLOBALFUNC::iserr(criterion[ii]))
            {
                int group_id = floor(criterion[ii]*group_total);
                if(group_id==group_total) group_id = group_total-1;
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
            ar & Alpha;


        }

        
        //globalallbegin
        //globalallend
        
        vector<vector<float>> Alpha;
        
        string member;
        vector<int> Ind;
        //constallbegin 
        const QL_MATRIX<QL_FLOAT> &close;
        const QL_CUBE<QL_FLOAT> &ask_high;
        const QL_CUBE<QL_FLOAT> &bid_high;
        const QL_CUBE<QL_FLOAT> &last_close;
        const QL_CUBE<QL_FLOAT> &last_mean;
        const QL_CUBE<QL_FLOAT> &mid_std;
        const QL_CUBE<QL_FLOAT> &mid_mean;
        const QL_CUBE<QL_FLOAT> &spread_high;
        const QL_CUBE<QL_FLOAT> &spread_low;
        const QL_CUBE<QL_FLOAT> &vol_sum;
        const QL_CUBE<QL_FLOAT> &ask_ave_price;
        const QL_CUBE<QL_FLOAT> &bid_ave_price;
        const QL_CUBE<QL_FLOAT> &sum_volume_price;
        const QL_MATRIX<QL_FLOAT> &buy_volume_small_order;
        const QL_MATRIX<QL_FLOAT> &close_net_inflow_rate_volume;
        //constallend
        
        bool flag_adj_update = false;

};
}

extern "C"
{
    AlphaBase * createStrategy(XMLCONFIG::Element *cfg)
    {
        AlphaBase * str = new alpha03(cfg);
        return str;
    }
}
//simplified