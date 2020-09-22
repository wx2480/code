#include "include/qlSim.h"
namespace
{
class dulin_alpha39 : public AlphaBase
{
    public:
    explicit dulin_alpha39(XMLCONFIG::Element *cfg):
        AlphaBase(cfg), 
        vwap(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.sum_volume_price")),
        close(dr.getData<QL_MATRIX<QL_FLOAT>>("adj_close")),
        close1(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.last_close")),
        low(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.last_low")),
        high(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.last_high"))
    {
        v_m.resize(GLOBAL::Dates.size());
        v_h.resize(GLOBAL::Dates.size());
        del1.resize(GLOBAL::Dates.size());
        for(int i = 0; i < GLOBAL::Dates.size(); ++i)
        {
            v_m[i].resize(GLOBAL::Instruments.size(),nan);
            v_h[i].resize(GLOBAL::Instruments.size(),nan);
            del1[i].resize(GLOBAL::Instruments.size(),nan);
        }
    }
    void generate(int di) override
    {
        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
        {
            if(valid[di][ii])
            {
                vector<float> c1(12,nan);
                vector<float> c2(12,nan);

                for(int i = 36; i < 48 ; ++i)
                {
                    c1[i - 36] = close1(di -delay, i ,ii);
                    c2[i - 36] = (low(di -delay,i,ii) - high(di -delay,i,ii));
                }
                v_m[di -delay][ii] = QL_Oputils::mean(c1);
                v_h[di -delay][ii] = QL_Oputils::mean(c2);
            }
        }
        QL_Oputils::rank(v_h[di -delay]);

         for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
        {
            if(valid[di][ii])
            {
                del1[di -delay][ii] = vwap(di -delay,48,ii) - vwap(di -delay,47,ii);
            }
        }
        vector<float> temp(GLOBAL::Instruments.size(),nan);
         for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
        {
            if(valid[di][ii])
            {
                temp[ii] = (close[di -delay][ii] - v_m[di -delay][ii]);
            }
        }
       QL_Oputils::rank(temp);

        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
        {
            if(valid[di][ii])
            {
                alpha[ii] =  -1*del1[di -delay][ii] * v_h[di -delay][ii] * temp[ii];
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
            ar & v_m;
            ar & v_h;
            ar & del1;
        }/*declare the data*/ 
        const QL_CUBE<QL_FLOAT> &vwap;
        const QL_MATRIX<QL_FLOAT> &close;
        const QL_CUBE<QL_FLOAT> &close1;
        const QL_CUBE<QL_FLOAT> &low;
        const QL_CUBE<QL_FLOAT> &high;

        vector<vector<float>> v_m;
        vector<vector<float>> v_h;
        vector<vector<float>> del1;

 
};
}
extern "C"
{
    AlphaBase * createStrategy(XMLCONFIG::Element *cfg)
    {
        AlphaBase * str = new dulin_alpha39(cfg);
        return str;
    }
}
