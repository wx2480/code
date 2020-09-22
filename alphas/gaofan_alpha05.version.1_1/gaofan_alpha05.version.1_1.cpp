
#include "include/qlSim.h"
namespace
{
class Gavin_Alpha_2nd_3 : public AlphaBase
{
    public:
    explicit Gavin_Alpha_2nd_3(XMLCONFIG::Element *cfg):
        AlphaBase(cfg), 
        volume(dr.getData<QL_MATRIX<QL_FLOAT>>("volume")), 
        vwap(dr.getData<QL_MATRIX<QL_FLOAT>>("adj_vwap")), 
        close(dr.getData<QL_MATRIX<QL_FLOAT>>("close")), 
        ret(dr.getData<QL_MATRIX<QL_FLOAT>>("adj_return")), 
        high(dr.getData<QL_MATRIX<QL_FLOAT>>("adj_high")), 
        Num1(cfg->getAttributeIntDefault("para1",5)),
        Num2(cfg->getAttributeIntDefault("para2",10))
    {
        MAVOL_5.resize(GLOBAL::Dates.size());
        MAVOL_10.resize(GLOBAL::Dates.size());
        MAVOL_5_10_delta.resize(GLOBAL::Dates.size());
        alpha_c.resize(GLOBAL::Dates.size());
        ret_c.resize(GLOBAL::Dates.size());
        eps.resize(GLOBAL::Dates.size());
        base_alpha.resize(GLOBAL::Dates.size());
        vwap_slop_5.resize(GLOBAL::Dates.size());
        vwap_slop_10.resize(GLOBAL::Dates.size());
        vwap_slop.resize(GLOBAL::Dates.size());
        high_slop.resize(GLOBAL::Dates.size());
        corr_vwap_volume.resize(GLOBAL::Dates.size());

        for(int di = 0; di < GLOBAL::Dates.size(); ++di)
        {
            MAVOL_5[di].resize(GLOBAL::Instruments.size(), nan);
            MAVOL_10[di].resize(GLOBAL::Instruments.size(), nan);
            MAVOL_5_10_delta[di].resize(GLOBAL::Instruments.size(), nan);
            alpha_c[di].resize(GLOBAL::Instruments.size(), nan);
            ret_c[di].resize(GLOBAL::Instruments.size(), nan);
            eps[di].resize(GLOBAL::Instruments.size(), nan);
            base_alpha[di].resize(GLOBAL::Instruments.size(), nan);
            vwap_slop_5[di].resize(GLOBAL::Instruments.size(), nan);
            vwap_slop_10[di].resize(GLOBAL::Instruments.size(), nan);
            vwap_slop[di].resize(GLOBAL::Instruments.size(), nan);
            high_slop[di].resize(GLOBAL::Instruments.size(), nan);
            corr_vwap_volume[di].resize(GLOBAL::Instruments.size(), nan);
        }

    }




    

    void generate(int di) override
    {
        QL_Oputils::qilin_ar_init(MAVOL_5, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(MAVOL_10, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(MAVOL_5_10_delta, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(alpha_c, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(ret_c, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(eps, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(base_alpha, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(vwap_slop_5, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(vwap_slop_10, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(vwap_slop, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(high_slop, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
        QL_Oputils::qilin_ar_init(corr_vwap_volume, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);

        calc_vwap_slop(di);
        calc_MAVOL(di);
        daily_calc(di);
        PID_control(di, 5, 0, 0, 1, base_alpha[di-delay]);

        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
        {
            if((valid[ di ][ ii ]))
            {
                alpha[ii] = base_alpha[di-delay][ii];
            }
        }
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



    void daily_calc(int di)
    {

        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
        {
            if((valid[ di ][ ii ]))
            {
                std::vector<float> past_MAVOL_5_10(2,nan);
                for(int i = 0; i < 2; ++ i)
                {
                    past_MAVOL_5_10[i] = MAVOL_5_10_delta[di-delay-i][ii]  ;

                }
                if(close[di-delay][ii] - vwap[di-delay][ii]>=0)
                {
                    base_alpha[di-delay][ii] = -ts_rank(past_MAVOL_5_10) ;
                }
                else
                {
                    base_alpha[di-delay][ii] = - (1 - ts_rank(past_MAVOL_5_10)) * vwap_slop[di-delay][ii] * 
                    (1 - corr_vwap_volume[di-delay][ii]) * high_slop[di-delay][ii] ;
                }
            }
        }
        return;
    }

    void calc_vwap_slop(int di)
    {
        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
        {
            std::vector<float> past_vwap_5(Num1,nan);
            std::vector<float> past_vwap_10(Num2,nan);
            std::vector<float> past_high_5(Num1,nan);
            std::vector<float> past_volume_10(Num2,nan);
            for(int i = 0; i < Num1; ++i)
            {
                past_vwap_5[i] = vwap[di-delay-i][ii];
                past_high_5[i] = high[di-delay-i][ii];
            }
            for(int i = 0; i < Num2; ++i)
            {
                past_vwap_10[i] = vwap[di-delay-i][ii];
                past_volume_10[i] = volume[di-delay-i][ii];
            }
            vwap_slop_5[di-delay][ii] = QL_Oputils::slop(past_vwap_5.begin(),past_vwap_5.end());
            vwap_slop_10[di-delay][ii] = QL_Oputils::slop(past_vwap_10.begin(),past_vwap_10.end());
            vwap_slop[di-delay][ii] = vwap_slop_5[di-delay][ii] - vwap_slop_10[di-delay][ii];
            high_slop[di-delay][ii] = QL_Oputils::slop(past_high_5.begin(),past_high_5.end());
            corr_vwap_volume[di-delay][ii] = QL_Oputils::corr(past_vwap_10,past_volume_10);
        }
        return;
    }

    void PID_control(int di, int cal_day, float Kp, float Ki, float Kd, vector<float>& base_alpha)
    {
    //Gloabal data: alpha_c, ret_c, eps;
        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ii) //loop for stocks set
        {
            ret_c[di-delay][ii] = ret[di-delay][ii];
        }
        QL_Oputils::rank(ret_c[di-delay]);
        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ii) //loop for stocks set
        {
            eps[di-delay][ii] = ret_c[di-delay][ii] - alpha_c[di-delay-1][ii];
            vector<float> past_eps(cal_day, nan);
            for (int i = 0; i < cal_day; ++i)
            {
                past_eps[i]  = eps[di-delay-i][ii];
            }

            float temp_kp = Kp*eps[di-delay][ii];
            if (GLOBALFUNC::iserr(temp_kp))
            {
                temp_kp = 0.0;
            }
            float temp_ki = Ki*QL_Oputils::decay_linear(past_eps);
            // float temp_ki = Ki*QL_Oputils::sum(past_eps);
             if (GLOBALFUNC::iserr(temp_ki))
            {
                temp_ki = 0.0;
            }
            float temp_kd = Kd*slop_ols(past_eps);
            if (GLOBALFUNC::iserr(temp_kd))
            {
                temp_kd = 0.0;
            }

            float deta_alpha = temp_kp + temp_ki + temp_kd;
            if (!GLOBALFUNC::iserr(deta_alpha))
            {
                base_alpha[ii] = base_alpha[ii]*(1.0 + deta_alpha);
            }
        }
        QL_Oputils::rank(base_alpha);

        for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ii)
        {
            alpha_c[di-delay][ii] = base_alpha[ii];
        }
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
            ar & MAVOL_5;
            ar & MAVOL_10;
            ar & MAVOL_5_10_delta;
            ar & alpha_c;
            ar & ret_c;
            ar & eps;
            ar & base_alpha;
            ar & vwap_slop_5;
            ar & vwap_slop_10;
            ar & vwap_slop;
            ar & high_slop;
            ar & corr_vwap_volume;

        }

        vector<vector<float>> MAVOL_5;
        vector<vector<float>> MAVOL_10;
        vector<vector<float>> MAVOL_5_10_delta;
        vector<vector<float>> alpha_c;
        vector<vector<float>> ret_c;
        vector<vector<float>> eps;
        vector<vector<float>> base_alpha;
        vector<vector<float>> vwap_slop_5;
        vector<vector<float>> vwap_slop_10;
        vector<vector<float>> vwap_slop;
        vector<vector<float>> high_slop;
        vector<vector<float>> corr_vwap_volume;


        const QL_MATRIX<QL_FLOAT> &volume;
        const QL_MATRIX<QL_FLOAT> &vwap;
        const QL_MATRIX<QL_FLOAT> &close;
        const QL_MATRIX<QL_FLOAT> &ret;
        const QL_MATRIX<QL_FLOAT> &high;
        int Num1;
        int Num2;
};
}
extern "C"
{
    AlphaBase * createStrategy(XMLCONFIG::Element *cfg)
    {
        AlphaBase * str = new Gavin_Alpha_2nd_3(cfg);
        return str;
    }
}

