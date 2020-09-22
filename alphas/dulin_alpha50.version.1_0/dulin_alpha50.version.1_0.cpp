#include "include/qlSim.h"
namespace
{
class dulin_alpha50 : public AlphaBase
{
    public:
    explicit dulin_alpha50(XMLCONFIG::Element *cfg):
        AlphaBase(cfg), 
        close1(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.last_close")),
        buy_large(dr.getData<QL_MATRIX<QL_FLOAT>>("wind.MoneyFlow.buy_value_large_order")),
        sell_large(dr.getData<QL_MATRIX<QL_FLOAT>>("wind.MoneyFlow.sell_value_large_order")),
        buy_medium(dr.getData<QL_MATRIX<QL_FLOAT>>("wind.MoneyFlow.buy_value_med_order")),
        sell_medium(dr.getData<QL_MATRIX<QL_FLOAT>>("wind.MoneyFlow.sell_value_med_order")),
        vwap(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.sum_volume_price")),
        buy_large_vol(dr.getData<QL_MATRIX<QL_FLOAT>>("wind.MoneyFlow.buy_volume_large_order")),
        sell_large_vol(dr.getData<QL_MATRIX<QL_FLOAT>>("wind.MoneyFlow.sell_volume_large_order")),
        buy_medium_vol(dr.getData<QL_MATRIX<QL_FLOAT>>("wind.MoneyFlow.buy_volume_med_order")),
        sell_medium_vol(dr.getData<QL_MATRIX<QL_FLOAT>>("wind.MoneyFlow.sell_volume_med_order")),
        Num1(cfg->getAttributeIntDefault("para1",30)),
        Num2(cfg->getAttributeIntDefault("para2",10))
    {
       
    }
    void generate(int di) override
    {
        vector<float> MA_10(GLOBAL::Instruments.size(),nan);
        vector<float> temp(GLOBAL::Instruments.size(),0);
        vector<float> MA_90(GLOBAL::Instruments.size(),nan);
        vector<float> temp1(GLOBAL::Instruments.size(),0);


        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
        {
            if(valid[di][ii])
            {
                vector<float> c(Num2,nan);
                for(int i = 0; i < Num2; ++i)
                {
                    c[i] = (buy_large[di -delay - i][ii] - sell_large[di - delay - i][ii])
                           + (buy_medium[di -delay - i][ii] - sell_medium[di - delay - i][ii]);
                }
                QL_Oputils::rank(c);
                MA_10[ii] = c[0];
            }
        }
        //QL_Oputils::rank(MA_10);
    
        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
        {
            if(valid[di][ii])
            {
                vector<float> c(Num1,nan);
                for(int i = 0; i < Num1; ++i)
                {
                    c[i] = (buy_large_vol[di -delay - i][ii] - sell_large_vol[di - delay - i][ii])
                        +  (buy_medium_vol[di -delay - i][ii] - sell_medium_vol[di - delay - i][ii]);
                }
                QL_Oputils::rank(c);
                MA_90[ii] = c[0];
            }
        }
       QL_Oputils::rank(MA_90);

        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
        {
            if(valid[di][ii])
            {
                //temp1[ii] = (vwap(di -delay,47,ii) + vwap(di -delay,48,ii));
                temp[ii] =  close1(di -delay,48,ii) - vwap(di - delay,42,ii);

                vector<float> c(24,nan);
                for(int i = 24; i < 49; ++i)
                {
                    c[i - 24] = vwap(di -delay,i,ii) - vwap(di -delay,i - 6,ii);
                }
                QL_Oputils::rank(c);
                temp1[ii] = c[23];
            }
        }
        QL_Oputils::rank(temp);
        QL_Oputils::rank(temp1);
 
        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
        {
            if(valid[di][ii])
            {  
                alpha[ii] = (1 - temp1[ii]) * (1 -temp[ii]) * (MA_10[ii]) * (MA_90[ii]);
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
        }/*declare the data*/ 
        const QL_CUBE<QL_FLOAT> &close1;
        const QL_MATRIX<QL_FLOAT> &buy_large;
        const QL_MATRIX<QL_FLOAT> &sell_large;
        const QL_MATRIX<QL_FLOAT> &buy_medium;
        const QL_MATRIX<QL_FLOAT> &sell_medium;
         const QL_CUBE<QL_FLOAT> &vwap;
         const QL_MATRIX<QL_FLOAT> &buy_large_vol;
         const QL_MATRIX<QL_FLOAT> &sell_large_vol;
         const QL_MATRIX<QL_FLOAT> &buy_medium_vol;
         const QL_MATRIX<QL_FLOAT> &sell_medium_vol;

        int Num1;
        int Num2;
 
};
}
extern "C"
{
    AlphaBase * createStrategy(XMLCONFIG::Element *cfg)
    {
        AlphaBase * str = new dulin_alpha50(cfg);
        return str;
    }
}
