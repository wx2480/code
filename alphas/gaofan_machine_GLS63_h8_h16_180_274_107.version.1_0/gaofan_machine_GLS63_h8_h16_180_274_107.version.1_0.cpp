
        #include "include/qlSim.h"
        namespace
        {
                        
class GLS63_h8_h16_180_274_107 : public AlphaBase
        {
        public:
                           
explicit GLS63_h8_h16_180_274_107(XMLCONFIG::Element *cfg):
                AlphaBase(cfg), 
                //ret(dr.getData<QL_MATRIX<QL_FLOAT>>("wind.MoneyFlow.moneyflow_pct_volume")), 
                ret(dr.getData<QL_MATRIX<QL_FLOAT>>("adj_vwap")), 
                //ret(dr.getData<QL_MATRIX<QL_FLOAT>>("volume")), 
                //ret(dr.getData<QL_MATRIX<QL_FLOAT>>("cap")),
                ret1(dr.getData<QL_MATRIX<QL_FLOAT>>("turnover")),
                //ret2(dr.getData<QL_VECTOR<QL_FLOAT>>("ZZ500.change")),
                ret3(dr.getData<QL_MATRIX<QL_FLOAT>>("turnover")),
relative_report_num_10d(dr.getData<QL_MATRIX<QL_FLOAT>> ("gogoal2.FocusStk.relative_report_num_10d")),

sum_volume_price(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.sum_volume_price")),

s_mfd_inflow_closevolume(dr.getData<QL_MATRIX<QL_FLOAT>>("wind.MoneyFlow.s_mfd_inflow_closevolume")),

adj_close(dr.getData<QL_MATRIX<QL_FLOAT>>("adj_close")),

ForecastRollRank_con_pb_rank(dr.getData<QL_MATRIX<QL_FLOAT>> ("gogoal2.ForecastRollRank.con_pb_rank")),

adj_open(dr.getData<QL_MATRIX<QL_FLOAT>>("adj_open")),

vwap(dr.getData<QL_MATRIX<QL_FLOAT>>("adj_vwap")),

volume_diff_large_trader(dr.getData<QL_MATRIX<QL_FLOAT>>("wind.MoneyFlow.volume_diff_large_trader")),

value_diff_small_trader(dr.getData<QL_MATRIX<QL_FLOAT>>("wind.MoneyFlow.value_diff_small_trader")),

sum_flow(dr.getData<QL_CUBE<QL_FLOAT>>("intervaldata_5min.sum_flow")),


                Num1(cfg->getAttributeIntDefault("para1",45)),
                Num2(cfg->getAttributeIntDefault("para2",45)),
                Num3(cfg->getAttributeIntDefault("para3",60)),
                Num4(cfg->getAttributeIntDefault("para4",60)),
                Num5(cfg->getAttributeIntDefault("para5",20))
                
                {
                     para1.resize(GLOBAL::Dates.size());
                for (int di=0;di<GLOBAL::Dates.size();++di)
                {
                    para1[di].resize ( GLOBAL::Instruments.size(), nan);
                }

                rank_alpha2.resize(GLOBAL::Dates.size());
                for (int di=0;di<GLOBAL::Dates.size();++di)
                {
                    rank_alpha2[di].resize ( GLOBAL::Instruments.size(), nan);
                }

                rank_alpha3.resize(GLOBAL::Dates.size());
                for (int di=0;di<GLOBAL::Dates.size();++di)
                {
                    rank_alpha3[di].resize ( GLOBAL::Instruments.size(), nan);
                }

                rank_alpha4.resize(GLOBAL::Dates.size());
                for (int di=0;di<GLOBAL::Dates.size();++di)
                {
                    rank_alpha4[di].resize ( GLOBAL::Instruments.size(), nan);
                }

                ret_yk.resize(GLOBAL::Dates.size());
                for (int di=0;di<GLOBAL::Dates.size();++di)
                {
                    ret_yk[di].resize ( GLOBAL::Instruments.size(), nan);
                }


                corr.resize(GLOBAL::Dates.size());
                for (int di=0;di<GLOBAL::Dates.size();++di)
                {
                    corr[di].resize ( GLOBAL::Instruments.size(), nan);
                }

                para2.resize(GLOBAL::Dates.size());
                for (int di=0;di<GLOBAL::Dates.size();++di)
                {
                    para2[di].resize ( GLOBAL::Instruments.size(), nan);
                }

                para3.resize(GLOBAL::Dates.size());
                for (int di=0;di<GLOBAL::Dates.size();++di)
                {
                    para3[di].resize ( GLOBAL::Instruments.size(), nan);
                }

                para4.resize(GLOBAL::Dates.size());
                for (int di=0;di<GLOBAL::Dates.size();++di)
                {
                    para4[di].resize ( GLOBAL::Instruments.size(), nan);
                }

                para5.resize(GLOBAL::Dates.size());
                for (int di=0;di<GLOBAL::Dates.size();++di)
                {
                    para5[di].resize ( GLOBAL::Instruments.size(), nan);
                }

                para6.resize(GLOBAL::Dates.size());
                for (int di=0;di<GLOBAL::Dates.size();++di)
                {
                    para6[di].resize ( GLOBAL::Instruments.size(), nan);
                }

                para7.resize(GLOBAL::Dates.size());
                for (int di=0;di<GLOBAL::Dates.size();++di)
                {
                    para7[di].resize ( GLOBAL::Instruments.size(), nan);
                }

                para8.resize(GLOBAL::Dates.size());
                for (int di=0;di<GLOBAL::Dates.size();++di)
                {
                    para8[di].resize ( GLOBAL::Instruments.size(), nan);
                }


                alpha_return_diff.resize(GLOBAL::Dates.size());
                for (int di=0;di<GLOBAL::Dates.size();++di)
                {
                    alpha_return_diff[di].resize ( GLOBAL::Instruments.size(), nan);
                }

                alpha_raw.resize(GLOBAL::Dates.size());
                for (int di=0;di<GLOBAL::Dates.size();++di)
                {
                    alpha_raw[di].resize ( GLOBAL::Instruments.size(), nan);
                }

                ret_c.resize(GLOBAL::Dates.size());
                for (int di=0;di<GLOBAL::Dates.size();++di)
                {
                    ret_c[di].resize ( GLOBAL::Instruments.size(), nan);
                }

                eps.resize(GLOBAL::Dates.size());
                for (int di=0;di<GLOBAL::Dates.size();++di)
                {
                    eps[di].resize ( GLOBAL::Instruments.size(), nan);
                }

                alpha_c.resize(GLOBAL::Dates.size());
                for (int di=0;di<GLOBAL::Dates.size();++di)
                {
                    alpha_c[di].resize ( GLOBAL::Instruments.size(), nan);
                }

                value_a.resize(GLOBAL::Dates.size());
                for (int di=0;di<GLOBAL::Dates.size();++di)
                {
                    value_a[di].resize ( GLOBAL::Instruments.size(), nan);
                }

                value_b.resize(GLOBAL::Dates.size());
                for (int di=0;di<GLOBAL::Dates.size();++di)
                {
                    value_b[di].resize ( GLOBAL::Instruments.size(), nan);
                }

                value_c.resize(GLOBAL::Dates.size());
                for (int di=0;di<GLOBAL::Dates.size();++di)
                {
                    value_c[di].resize ( GLOBAL::Instruments.size(), nan);
                }

               base_alpha.resize(GLOBAL::Dates.size());
               for(int di = 0; di < GLOBAL::Dates.size();++di)
               {
                    base_alpha[di].resize(GLOBAL::Instruments.size(),nan);
               }

               binary_state.resize(GLOBAL::Dates.size());
               for(int di = 0; di < GLOBAL::Dates.size();++di)
               {
                    binary_state[di].resize(GLOBAL::Instruments.size(),nan);
               }


                   

                }
                
            void generate(int di) override
            {

                QL_Oputils::qilin_ar_init(para1, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
                QL_Oputils::qilin_ar_init(rank_alpha2, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
                QL_Oputils::qilin_ar_init(rank_alpha3, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
                QL_Oputils::qilin_ar_init(rank_alpha4, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
                QL_Oputils::qilin_ar_init(ret_yk, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
                QL_Oputils::qilin_ar_init(corr, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
                QL_Oputils::qilin_ar_init(para2, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);

                QL_Oputils::qilin_ar_init(para3, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
                QL_Oputils::qilin_ar_init(para4, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
                QL_Oputils::qilin_ar_init(para5, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
                QL_Oputils::qilin_ar_init(para6, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
                QL_Oputils::qilin_ar_init(para7, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
                QL_Oputils::qilin_ar_init(para8, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
                QL_Oputils::qilin_ar_init(alpha_return_diff, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
                QL_Oputils::qilin_ar_init(alpha_raw, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);

                QL_Oputils::qilin_ar_init(ret_c, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
                QL_Oputils::qilin_ar_init(eps, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
                QL_Oputils::qilin_ar_init(alpha_c, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);


                QL_Oputils::qilin_ar_init(value_a, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
                QL_Oputils::qilin_ar_init(value_b, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
                QL_Oputils::qilin_ar_init(value_c, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);

                QL_Oputils::qilin_ar_init(base_alpha, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);
                QL_Oputils::qilin_ar_init(binary_state, GLOBAL::Dates.size(), GLOBAL::Instruments.size(), nan);

                daily_calculate(di);

                

                for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
                {   
                    if((valid[ di ][ ii ]))
                    {
                        alpha[ii]=rank_alpha2[di-delay][ii];
                        //alpha[ii]=pow((1.0+rank_alpha2[di-delay][ii]),(1.0+para5[di-delay][ii]));

                    }
                }

                float mean_alpha = QL_Oputils::median(alpha);
                for(int ii = 0; ii < GLOBAL::Instruments.size();++ii)
                {
                    if(valid[di][ii])
                    {
                        if(GLOBALFUNC::iserr(alpha[ii]))
                        {
                            alpha[ii] = mean_alpha;
                        }
                    }
                }


                 return;
            }

            void daily_calculate(int di)
            {
                for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
                {   
                    if((valid[ di ][ ii ]))
                    {
                    /*para6[di-delay][ii]=ret[di-delay][ii];
                    vector<float> H70_cor_gestock_corr(10,nan);
                    vector<float> H70_cor_dapanindex_corr(10,nan);
                    for(int i=0;i<10;i++)
                    {
                    H70_cor_gestock_corr[i]=ret3[di-delay-i][ii];
                    H70_cor_dapanindex_corr[i]=ret2[di-delay-i];
                    }
                    para7[di-delay][ii]=(QL_Oputils::corr(H70_cor_gestock_corr,H70_cor_dapanindex_corr));*/
                    
                           
vector<float> price_diversion(10,nan);
                 for(int i = 0; i < 10; ++i)
                 {
                   price_diversion[i] = adj_close(di-delay-i,ii) - vwap(di-delay-i,ii);
                 }
                  para1[di-delay][ii]=-ts_rank(price_diversion);

vector<float> corr_sum_volume_price_43_48(6,nan);
                vector<float> corr_sum_flow_43_48_sum_volume_price(6,nan);
                for(int i=0;i<6;i++)
                {
                    corr_sum_volume_price_43_48[i]=sum_volume_price(di-delay,48-i,ii);
                    corr_sum_flow_43_48_sum_volume_price[i]=sum_flow(di-delay,48-i,ii);
                }
                para2[di-delay][ii]=-QL_Oputils::corr(corr_sum_volume_price_43_48,corr_sum_flow_43_48_sum_volume_price);//0.7,6.11

vector<float> L180_cnt_tsrankmoneyflow_5(5,nan);
 for(int i=0;i<5;i++)
{
L180_cnt_tsrankmoneyflow_5[i]=volume_diff_large_trader[di-delay-i][ii];
}
para3[di-delay][ii]=(1-ts_rank(L180_cnt_tsrankmoneyflow_5));


vector<float> L274_cnt_open_corr(10,nan);
vector<float> L274_cnt_pcorml_corr(10,nan);
for(int i=0;i<10;i++)
{
L274_cnt_open_corr[i]=adj_open[di-delay-i][ii];
L274_cnt_pcorml_corr[i]=value_diff_small_trader[di-delay-i-1][ii];
}
para4[di-delay][ii]=-(QL_Oputils::corr(L274_cnt_open_corr,L274_cnt_pcorml_corr));


std::vector<float> L107_past_report_37(20,nan);
            std::vector<float> L107_past_eps_37(20,nan);
            for(int i = 0; i < 20 ; i++)
            {
                L107_past_report_37[i] = relative_report_num_10d[di-delay-i][ii];
                L107_past_eps_37[i] = ForecastRollRank_con_pb_rank[di-delay-i][ii];
            }
            para5[di-delay][ii] = QL_Oputils::corr(L107_past_report_37,L107_past_eps_37); //0.2-1 0.1-117 0.24
            


                    }

                }
                
                QL_Oputils::rank(para1[di-delay]);
                QL_Oputils::rank(para2[di-delay]);
                QL_Oputils::rank(para3[di-delay]);
                QL_Oputils::rank(para4[di-delay]);
                QL_Oputils::rank(para5[di-delay]);
                //QL_Oputils::rank(para6[di-delay]);
                //QL_Oputils::rank(para7[di-delay]);
                
                for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
                {   
                    if((valid[ di ][ ii ]))
                    {
                        
                        rank_alpha3[di-delay][ii]=pow(1.0+para1[di-delay][ii],1.0+(Num1/100.0))*pow(1.0+para2[di-delay][ii],1.0+(Num2/100.0))*pow(1.0+para3[di-delay][ii],1.0+(Num3/100.0))*pow(1.0+para4[di-delay][ii],1.0+(Num4/100.0));
                        //rank_alpha3[di-delay][ii]=pow((1.0+para1[di-delay][ii]),para3[di-delay][ii])*pow((1.0+para2[di-delay][ii]),para4[di-delay][ii]);
                        //rank_alpha3[di-delay][ii]=pow(1.0+para1[di-delay][ii],1.0+(Num1/100.0))*pow(1.0+para2[di-delay][ii],1.0+(Num2/100.0));
                    }
                }
                
                weighted_group_rank(di, 100, rank_alpha3[di-delay]);

                QL_Oputils::rank(rank_alpha3[di-delay]);

                for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ ii)
                {   
                    if((valid[ di ][ ii ]))
                    {
                        vector<float> rank_alpha_10(Num5,nan);
                        for(int i=0;i<Num5;i++)
                        {
                            if(!GLOBALFUNC::iserr(para5[di-delay-i][ii]))
                            {
                            rank_alpha_10[i]=pow((1.0+rank_alpha3[di-delay-i][ii]),(1.0+para5[di-delay-i][ii]));
                            }
                            else
                            {
                            rank_alpha_10[i]=pow((1.0+rank_alpha3[di-delay-i][ii]),(1.0+0.0));
                            }
                            //rank_alpha_10[i]=rank_alpha3[di-delay-i][ii];

                        }
                       
                       rank_alpha2[di-delay][ii]=ts_rank(rank_alpha_10);

                    }

                }

                //QL_Oputils::rank(rank_alpha2[di-delay]);

                return;
            }

            
            void weighted_group_rank(int di, int group_total, vector<float>&cross_sectional_vect) 
            {

                vector<float> criterion(GLOBAL::Instruments.size(),nan);

                for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ii)
                {
                    
criterion[ii] =s_mfd_inflow_closevolume[di-delay][ii];


                }
                QL_Oputils::rank(criterion);

                
                
                vector<vector<float>> grouped_data(group_total, vector<float>(GLOBAL::Instruments.size(), nan));
                
                vector<vector<float>> weight_standards(group_total, vector<float>(GLOBAL::Instruments.size(), nan));

                vector<vector<float>> weight_standards1(group_total, vector<float>(GLOBAL::Instruments.size(), nan));

                for(int ii = 0; ii < GLOBAL::Instruments.size(); ++ii)
                {
                    if(!GLOBALFUNC::iserr(criterion[ii]))
                    {
                        int group_id = floor(criterion[ii] * group_total);

                        if(group_id == group_total) group_id = group_total-1;

                        grouped_data[group_id][ii] = cross_sectional_vect[ii];

                        //vector<float> price_diversion_w(10,nan);
                        //for(int i = 0; i < 10; ++i)
                        //{
                        //price_diversion_w[i] = ret1[di-delay-i][ii] - ret[di-delay-i][ii];
                        //}
                        //weight_standards[group_id][ii]=-ts_rank(price_diversion_w);

                        

                        weight_standards[group_id][ii]= ret[di-delay][ii];

                        /*vector<float> hhb_cnt_vol_slop(10,nan);
                        for(int i = 0; i < 10; ++ i){hhb_cnt_vol_slop[i] = ret[di-delay-i][ii];}

                        weight_standards[group_id][ii]=(QL_Oputils::decay_exp(hhb_cnt_vol_slop,0.25));*/
                        
                        
                        

                        

                     }
                 }
                 vector<float> group_weight(group_total, nan);

                 for(int group_id = 0; group_id < group_total; ++ group_id)
                 {
                    QL_Oputils::rank(grouped_data[group_id]);

                    //group_weight[group_id] = QL_Oputils::mean(weight_standards[group_id]);

                    group_weight[group_id] = QL_Oputils::median(weight_standards[group_id]);

                    //group_weight[group_id] = QL_Oputils::hhv(weight_standards[group_id].begin(),weight_standards[group_id].end());

                    //group_weight[group_id] = QL_Oputils::llv(weight_standards[group_id].begin(),weight_standards[group_id].end());

                    //group_weight[group_id] = (QL_Oputils::median(weight_standards[group_id])+QL_Oputils::hhv(weight_standards[group_id].begin(),weight_standards[group_id].end()))/2.0;
                    
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


            
            
            float ts_rank(vector<float> &vect) 
            {
                QL_Oputils::rank(vect);
                return vect[0];
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
                    /*initialize global parameters here*/
                    ar & para1;
                    ar & rank_alpha2;
                    ar & rank_alpha3;
                    ar & rank_alpha4;
                    ar & ret_yk;
                    ar & corr;
                    ar & para2;
                    ar & para3;
                    ar & para4;
                    ar & para5;
                    ar & para6;
                    ar & para7;
                    ar & para8;
                    ar & alpha_return_diff;
                    ar & alpha_raw;

                    ar & ret_c;
                    ar & eps;
                    ar & alpha_c;

                    ar & value_a;
                    ar & value_b;
                    ar & value_c;

                    ar & base_alpha;
                    ar & binary_state;
           
                    
                }
                 
                vector<vector<float> > para1;
                vector<vector<float> > rank_alpha2;
                vector<vector<float> > rank_alpha3;
                vector<vector<float> > rank_alpha4;
                vector<vector<float> > ret_yk;
                vector<vector<float> > corr;
                vector<vector<float> > para2;
                vector<vector<float> > para3;
                vector<vector<float> > para4;
                vector<vector<float> > para5;
                vector<vector<float> > para6;
                vector<vector<float> > para7;
                vector<vector<float> > para8;
                vector<vector<float> > alpha_return_diff;
                vector<vector<float> > alpha_raw;

                vector<vector<float> > ret_c;
                vector<vector<float> > eps;
                vector<vector<float> > alpha_c;

                vector<vector<float>> value_a;
                vector<vector<float>> value_b;
                vector<vector<float>> value_c;

                vector<vector<float>> base_alpha;
                vector<vector<float>> binary_state;
                
                
                const QL_MATRIX<QL_FLOAT> &ret;
                const QL_MATRIX<QL_FLOAT> &ret1;
                //const QL_VECTOR<QL_FLOAT> &ret2;
                const QL_MATRIX<QL_FLOAT> &ret3;
const QL_MATRIX<QL_FLOAT> &relative_report_num_10d;

const QL_CUBE<QL_FLOAT> &sum_volume_price;

const QL_MATRIX<QL_FLOAT> &s_mfd_inflow_closevolume;

const QL_MATRIX<QL_FLOAT> &adj_close;

const QL_MATRIX<QL_FLOAT> &ForecastRollRank_con_pb_rank;

const QL_MATRIX<QL_FLOAT> &adj_open;

const QL_MATRIX<QL_FLOAT> &vwap;

const QL_MATRIX<QL_FLOAT> &volume_diff_large_trader;

const QL_MATRIX<QL_FLOAT> &value_diff_small_trader;

const QL_CUBE<QL_FLOAT> &sum_flow;

 
                int Num1;
                int Num2;
                int Num3;
                int Num4;
                int Num5;
            
               
        };
        }
        extern "C"
        {
            AlphaBase * createStrategy(XMLCONFIG::Element *cfg)
            {

                           
AlphaBase * str = new GLS63_h8_h16_180_274_107(cfg);
                        return str;
            }
        }
                           
