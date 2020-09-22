#include "include/qlSim.h"
namespace
{
class alpha40 : public AlphaBase
{
    public:
    explicit alpha40(XMLCONFIG::Element *cfg):
        AlphaBase(cfg), 
        close(dr.getData<QL_MATRIX<QL_FLOAT>>("adj_close")),
        amount(dr.getData<QL_MATRIX<QL_FLOAT>>("volume")),
        vwap(dr.getData<QL_MATRIX<QL_FLOAT>>("vwap"))
    {
        ret.resize(GLOBAL::Dates.size());
  
        for(int i = 0; i < GLOBAL::Dates.size(); ++i)
        {
            ret[i].resize(GLOBAL::Instruments.size(),nan);
        }
    }
    void generate(int di) override
    {
        vector<float> liqu(GLOBAL::Instruments.size(),nan);
        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
        {
            if(valid[di][ii])
            {
                if(abs( log(close[di -delay][ii]) - log(close[di -delay - 1][ii])) > 0)
                {
                    liqu[ii] = amount[di -delay][ii] / ( log(close[di -delay][ii]) - log(close[di -delay - 1][ii]) );
                }
                ret[di -delay][ii] =  ( log(close[di -delay][ii]) - log(close[di -delay - 1][ii]) );
            }
        }
        QL_Oputils::rank(liqu);

        int flag1 = 0;
        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
        {
            if(valid[di][ii])
            {
                if(liqu[ii] > 0.7)
                {
                    flag1++;
                }
            }
        }

        vector<float> temp(flag1,nan);
        int flag2 = 0;

        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
        {
            if(valid[di][ii])
            {
                if(liqu[ii] > 0.7)
                {
                    temp[flag2] = log(close[di -delay][ii]) - log(close[di -delay -1][ii]);
                    flag2++;
                }
            }
        }

        float kur = QL_Oputils::kurtosis(temp);
        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
        {
            if(valid[di][ii])
            {                    
                if(kur > 3)
                {
                    alpha[ii] = liqu[ii];
                }
                else
                {
                    alpha[ii] = -1*ret[di -delay][ii] * liqu[ii];
                }
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
            ar & ret;
        }
        const QL_MATRIX<QL_FLOAT> &close;   /*declare the data*/
        const QL_MATRIX<QL_FLOAT> &amount;
        const QL_MATRIX<QL_FLOAT> &vwap;

        vector<vector<float>> ret;
     
};
}
extern "C"
{
    AlphaBase * createStrategy(XMLCONFIG::Element *cfg)
    {
        AlphaBase * str = new alpha40(cfg);
        return str;
    }
}
