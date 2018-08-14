set search_path to mimiciii;

-- step 1. Find adults between ages 16 and 89.
-- step 2. Create chartevents table with only adults.
-- step 3. Update table to include patients from icu. This includes MICU, SICU, CCU, CSRU.
-- step 4. Query chartevents_adults for patients that have expierenced a hypotensive episode. 
-- step 5. Remove Patients that are coded with comfort measures only (CMO).
-- step 6. Find Patients pre-serum creatinine and post-serum creatinine values.
-- step 7. Obtain Patients that recieved Fluids treatment.
-- step 8. Repeat step 7 for patients given vasopressors.
-- step 9. Consolidate steps 7 and 8 into 1 table.

/* Remove patients younger than 16. */
WITH first_admission_time AS (
SELECT 
    p.subject_id, p.dob, p.gender, 
    MIN (a.admittime) AS first_admittime
FROM patients p
INNER JOIN admissions a
ON p.subject_id = a.subject_id
GROUP BY p.subject_id, p.dob, p.gender, a.hadm_id
ORDER BY a.hadm_id, p.subject_id
),
age AS (
SELECT 
    subject_id,  dob, gender, first_admittime, 
    age(first_admittime, dob) 
        AS first_admit_age, 
    CASE
       WHEN age(first_admittime, dob) > '89 years'
            THEN '>89'
        WHEN age(first_admittime, dob) >= '15 years'
            THEN 'adult'
        WHEN age(first_admittime, dob) <= '1  year'
            THEN 'neonate'
        ELSE 'middle'
        END AS age_group
FROM first_admission_time
ORDER BY subject_id
)
SELECT *
into subject_age
FROM age;

-- step 2 and step 3. Find adults in the ICU's: MICU, CCU, CSRU, SICU.
select *
into chartevents_adult
from chartevents
where chartevents.subject_id in (
  select subject_age.subject_id
  from subject_age
  where age_group = 'adult'
	) 
  and chartevents.icustay_id in (
  select icustay_id
  from icustays
  where first_careunit in ('MICU','CCU','SICU','CSRU') or last_careunit in ('MICU','CCU','SICU','CSRU')  
);

-- step 4.
/* Extract entire Hypotensive cohort from Chartevents table */

/* The beginning of an HE was deﬁned as the time of the ﬁrst of two consecutive MAP measurements 
<60 mm Hg, preceded by two consecutive MAP values above 60 mm Hg. The end of an HE was deﬁned as 
the time of the ﬁrst of two consecutive MAP measurements >60 mm Hg, preceded by two consecutive MAP 
values #60 mm Hg. A MAP threshold of 60 mm Hg was used in the study based on the ﬁnding that autoregulation 
ceases and blood ﬂow becomes pressure dependent below this level in various organs. */


/* MAP(ABPm) in metavision - 220052
   MAP in Carevue - 52,
   IABP mean in Carevue - 224,
   Manual BP mean in Carevue - 443,
   Arterial BP mean #2 in Carevue - 6702,
   Arterial Mean #3 in Carevue - 6927 */
  
 
with HE_CE_adult_START AS (
select chartevents_adult.icustay_id, itemid, valuenum, next_charttime as he_start
from (select chartevents_adult.icustay_id, chartevents_adult.charttime, chartevents_adult.itemid, chartevents_adult.valuenum,
             row_number()  over (partition by chartevents_adult.icustay_id order by chartevents_adult.charttime) as rn,
             lag(chartevents_adult.valuenum,2) over (partition by chartevents_adult.icustay_id order by chartevents_adult.charttime) as prev_value2,
             lag(chartevents_adult.valuenum) over (partition by chartevents_adult.icustay_id order by chartevents_adult.charttime) as prev_value,
             lead(chartevents_adult.valuenum) over (partition by chartevents_adult.icustay_id order by chartevents_adult.charttime) as next_value,
      		 lead(chartevents_adult.charttime) over (partition by chartevents_adult.icustay_id order by chartevents_adult.charttime) as next_charttime
      from chartevents_adult
      where chartevents_adult.itemid in (225312, 220052, 6926, 52, 6702, 456, 220181)
     ) chartevents_adult
where chartevents_adult.valuenum <= 60 and next_value < 60 and 
      prev_value > 60 and prev_value2 > 60 and (next_charttime - chartevents_adult.charttime) <> '00:00:00'
) 
 , HE_CE_adult_END AS(
        select chartevents_adult.icustay_id, (last_charttime - he_start) as HE_duration, he_start, last_charttime as HE_end
        from (select chartevents_adult.icustay_id, chartevents_adult.charttime, chartevents_adult.itemid, chartevents_adult.valuenum,
                     row_number()  over (partition by chartevents_adult.icustay_id order by chartevents_adult.charttime) as rn2,
                     lag(chartevents_adult.valuenum,2) over (partition by chartevents_adult.icustay_id order by chartevents_adult.charttime) as prev_value_new2,
                     lag(chartevents_adult.valuenum) over (partition by chartevents_adult.icustay_id order by chartevents_adult.charttime) as prev_value_new,
                     lead(chartevents_adult.valuenum) over (partition by chartevents_adult.icustay_id order by chartevents_adult.charttime) as final_value,
                     lead(chartevents_adult.charttime) over (partition by chartevents_adult.icustay_id order by chartevents_adult.charttime) as last_charttime
              from chartevents_adult
              where chartevents_adult.itemid in (225312, 220052, 6926, 52, 6702, 456, 220181)
             ) chartevents_adult, HE_CE_adult_START
        where (chartevents_adult.charttime - HE_CE_adult_START.he_start) between '00:00:00' and '24:00:00' and chartevents_adult.valuenum > 60 and final_value > 60 and
              prev_value_new <= 60 and prev_value_new2 <= 60 and (last_charttime - chartevents_adult.charttime) <> '00:00:00'
) 
 , HE_Cohort AS(
        select distinct(HE_CE_adult_END.icustay_id), HE_duration, HE_CE_adult_END.he_start, HE_CE_adult_END.he_end 
        from HE_CE_adult_END
        inner join HE_CE_adult_START on HE_CE_adult_END.icustay_id = HE_CE_adult_START.icustay_id
)
select HE_Cohort.icustay_id, min(HE_Cohort.HE_duration) as he_duration, min(HE_Cohort.he_start) as he_starttime, min(HE_Cohort.he_end) as he_endtime, intime as admittime
into HE_Cohort
from HE_Cohort, icustays
where (HE_Cohort.he_end - icustays.intime) between '00:00:00' and '48:00:00'
    and HE_Cohort.icustay_id = icustays.icustay_id 
group by HE_Cohort.icustay_id, intime
order by HE_Cohort.icustay_id  
;

-- This code finds the first HE epsiode within first 48 hours of admission.
-- Values for MAP to get HE episode.

/* MAP(ABPm) in metavision - 220052
   MAP in Carevue - 52,
   IABP mean in Carevue - 224,
   Manual BP mean in Carevue - 443,
   Arterial BP mean #2 in Carevue - 6702,
   Arterial Mean #3 in Carevue - 6927 */
  
-- step 5. Remove pateints with CMO's.
WITH icustay_CMO AS (
    select HE_Cohort.icustay_id, chartevents.itemid, chartevents.value, chartevents.charttime, HE_Cohort.he_duration, HE_Cohort.he_starttime, HE_Cohort.he_endtime, 
        HE_Cohort.admittime 
         from chartevents
         inner join HE_Cohort on chartevents.icustay_id = HE_Cohort.icustay_id
    		)
	, icustay_CMO_id AS (
        select *
            from icustay_CMO
            where  value in ('Comfort Measures','Comfort measures only','CMO') and
                        (he_starttime - charttime ) between '-24:00:00' and '24:00:00' or
                   value in ('Comfort Measures','Comfort measures only','CMO') and
                        (he_endtime - charttime ) between '-24:00:00' and '24:00:00'
        ) 
select distinct on (1) HE_Cohort.icustay_id, HE_Cohort.he_duration, HE_Cohort.he_starttime, HE_Cohort.he_endtime, HE_Cohort.admittime 
into HE_Cohort_final
from HE_Cohort, icustay_CMO_id
where HE_Cohort.icustay_id not in (select icustay_id from icustay_CMO_id);

drop table HE_Cohort;

-- step 6 Get response variable. 
-- This code finds the serum creatinine values 24hr before HE onset.
WITH CE_HE_Small AS (
    select HE_Cohort_final.icustay_id, chartevents.itemid, chartevents.valuenum, chartevents.charttime, HE_Cohort_final.he_duration, HE_Cohort_final.he_starttime, HE_Cohort_final.he_endtime, 
        HE_Cohort_final.admittime 
         from chartevents
         inner join HE_Cohort_final on chartevents.icustay_id = HE_Cohort_final.icustay_id
    		)
	, Pre_Serum_CE AS (
    select icustay_id, itemid, valuenum, charttime, he_starttime, he_endtime, admittime 
    from CE_HE_Small
    where itemid in (791, 1525, 220615, 3750) and 
          ( he_starttime - charttime ) between '00:00:00' and '24:00:00'
    order by icustay_id,itemid, valuenum, he_starttime, admittime
        )
    , CE_HE_PS AS (    
    SELECT distinct PS.icustay_id, PS.valuenum as pre_se_cr, he_starttime, he_endtime, admittime 
    FROM Pre_Serum_CE PS
    INNER JOIN
        (SELECT icustay_id, MAX(charttime) AS MaxChartTime
        FROM Pre_Serum_CE
        GROUP BY icustay_id) groupedPS 
    ON PS.icustay_id = groupedPS.icustay_id 
    AND PS.charttime = groupedPS.MaxChartTime
	)
select HE_Cohort_final.icustay_id, HE_Cohort_final.he_starttime, HE_Cohort_final.he_endtime, HE_Cohort_final.admittime, pre_se_cr
into pre_sc_HE
from HE_Cohort_final
LEFT JOIN CE_HE_PS
ON HE_Cohort_final.icustay_id = CE_HE_PS.icustay_id;

-- This code finds the serum creatinine values 72hr after HE onset.
WITH CE_HE_Small AS (
    select HE_Cohort_final.icustay_id, chartevents.itemid, chartevents.valuenum, chartevents.charttime, HE_Cohort_final.he_duration, HE_Cohort_final.he_starttime, HE_Cohort_final.he_endtime, 
        HE_Cohort_final.admittime 
         from chartevents
         inner join HE_Cohort_final on chartevents.icustay_id = HE_Cohort_final.icustay_id
    		)
	, Post_Serum_CE AS (
    select icustay_id, itemid, valuenum, charttime, he_starttime, he_endtime, admittime 
    from CE_HE_Small
    where itemid in (791, 1525, 220615, 3750) and 
          ( charttime - he_endtime ) between '00:00:00' and '72:00:00'
    order by icustay_id,itemid, valuenum, he_starttime, admittime
        )
    , CE_HE_PS AS (    
    SELECT distinct PS.icustay_id, PS.valuenum as post_se_cr, he_starttime, he_endtime, admittime 
    FROM Post_Serum_CE PS
    INNER JOIN
        (SELECT icustay_id, MAX(charttime) AS MaxChartTime
        FROM Post_Serum_CE
        GROUP BY icustay_id) groupedPS 
    ON PS.icustay_id = groupedPS.icustay_id 
    AND PS.charttime = groupedPS.MaxChartTime
	)
select pre_sc_HE.icustay_id, pre_sc_HE.he_starttime, pre_sc_HE.he_endtime, pre_sc_HE.admittime, pre_se_cr, post_se_cr
into post_se_HE
from pre_sc_HE
LEFT JOIN CE_HE_PS
ON pre_sc_HE.icustay_id = CE_HE_PS.icustay_id;

-- creates response value: change in serum creatinine. 
select *, (he_endtime - he_starttime) as he_duration, (post_se_cr - pre_se_cr) as chg_ser
into HE_ds1
from post_se_HE;

drop table post_se_HE;
drop table pre_sc_HE;

--to find itemid for serum creatinine.
select * from d_items
WHERE label like '%serum%'
OR label like '%creatinine%';


-- step 7. Patients with Fluids Treatment.

--patients in cv.
WITH fluid_patients_cv AS (
     select he_ds1.icustay_id , inputevents_cv.itemid as cvfluid_id, inputevents_cv.charttime as cvfluid_input_time, inputevents_cv.amount as cvfluid_vol
     from he_ds1
     inner join inputevents_cv on he_ds1.icustay_id = inputevents_cv.icustay_id
     where inputevents_cv.charttime BETWEEN he_ds1.he_starttime AND he_ds1.he_endtime
            AND ( inputevents_cv.itemid in
				(
                    30176, 30161, 30020, 30018, 30353, 30352, 30143, 45232, 30160, 30014, 30015, 30060, 225825, 45374, 
                    45358, 30061, 44399, 221001, 221003, 220966, 220963, 220965, 220964, 226401, 220979, 221195,45532, 
                    227344,1634, 30021,5455, 6190,4647, 44440,44053, 43354,30168, 1494, 221212, 220958, 220959, 220960, 
                    220961, 220954, 220962, 221213, 5333
					) )
        		AND amount > 250
        ) 
select distinct on (1) he_ds1.*, fluid_patients_cv.cvfluid_id, fluid_patients_cv.cvfluid_input_time, fluid_patients_cv.cvfluid_vol 
-- multiple hits within HE episode so only count it as one. 
into he_ds_fluidscv
from he_ds1
LEFT JOIN fluid_patients_cv
ON  he_ds1.icustay_id = fluid_patients_cv.icustay_id; 

-- patients in mv.    
WITH fluid_patients_mv AS (
        select he_ds1.* , inputevents_mv.itemid as mvfluid_id, inputevents_mv.starttime as mvfluid_input_time, inputevents_mv.amount as mvfluid_vol ,
    		   inputevents_mv.cancelreason, inputevents_mv.icustay_id as icd
            from he_ds1
        	inner join inputevents_mv on he_ds1.icustay_id = inputevents_mv.icustay_id
        		where inputevents_mv.starttime BETWEEN he_ds1.he_starttime AND he_ds1.he_endtime
                AND ( inputevents_mv.itemid in
				(
                    30176, 30161, 30020, 30018, 30353, 30352, 30143, 45232, 30160, 30014, 30015, 30060, 225825, 45374, 
                    45358, 30061, 44399, 221001, 221003, 220966, 220963, 220965, 220964, 226401, 220979, 221195,45532, 
                    227344,1634, 30021,5455, 6190,4647, 44440,44053, 43354,30168, 1494, 221212, 220958, 220959, 220960, 
                    220961, 220954, 220962, 221213, 5333
					) )
        		AND amount > 250
    			AND cancelreason = 0
        ) 
select distinct on (1) he_ds_fluidscv.*, fluid_patients_mv.mvfluid_id, fluid_patients_mv.mvfluid_input_time, fluid_patients_mv.mvfluid_vol 
into he_ds_fluids_all
from he_ds_fluidscv
LEFT JOIN  fluid_patients_mv
ON he_ds_fluidscv.icustay_id = fluid_patients_mv.icustay_id; 

-- view variables to figure out what to keep.
select * from he_ds_fluids_all limit 10;

-- now make indicator variable if patient recieved fluids treatment.
select *, ( CASE
		WHEN cvfluid_vol is not Null or mvfluid_vol is not Null Then 1
    	ELSE 0
    	END ) AS IV_fluid_ind
into HE_ds_fluids_trt
from he_ds_fluids_all;

drop table he_ds_fluids_all;
drop table  he_ds_fluidscv;

-- step 8. Repeat step 7 for vasopressors.

--patients in cv.
WITH vp_patients_cv AS (
     select he_ds1.icustay_id , inputevents_cv.itemid as cv_vp_id, inputevents_cv.charttime as cvvp_input_time, inputevents_cv.amount as cvvp_amount
     from he_ds1
     inner join inputevents_cv on he_ds1.icustay_id = inputevents_cv.icustay_id
     where inputevents_cv.charttime BETWEEN he_ds1.he_starttime AND he_ds1.he_endtime
            AND ( inputevents_cv.itemid in
				(
                    30047,30120,221906, 30044,30119,30309,221289,30127,30128,221749,
                  30051,222315, 30043,30307,221662, 30042,30306,221653, 30125,221986 
					) )
        ) 
select distinct on (1) he_ds1.*, vp_patients_cv.cv_vp_id, vp_patients_cv.cvvp_input_time, vp_patients_cv.cvvp_amount
-- multiple hits within HE episode so only count it as one. 
into he_ds_vpcv
from he_ds1
LEFT JOIN vp_patients_cv
ON  he_ds1.icustay_id = vp_patients_cv.icustay_id; 

-- patients in mv.    
WITH vp_patients_mv AS (
        select he_ds1.* , inputevents_mv.itemid as mvvp_id, inputevents_mv.starttime as mvvp_input_time, inputevents_mv.amount as mvvp_amount ,
    		   inputevents_mv.cancelreason, inputevents_mv.icustay_id as icd
            from he_ds1
        	inner join inputevents_mv on he_ds1.icustay_id = inputevents_mv.icustay_id
        		where inputevents_mv.starttime BETWEEN he_ds1.he_starttime AND he_ds1.he_endtime
                AND ( inputevents_mv.itemid in
				(
                    30047,30120,221906, 30044,30119,30309,221289,30127,30128,221749,
                  30051,222315, 30043,30307,221662, 30042,30306,221653, 30125,221986 
					) )
    			AND cancelreason = 0
        ) 
select distinct on (1) he_ds_vpcv.*, vp_patients_mv.mvvp_id, vp_patients_mv.mvvp_input_time, vp_patients_mv.mvvp_amount 
into he_ds_vp_all
from he_ds_vpcv
LEFT JOIN  vp_patients_mv
ON he_ds_vpcv.icustay_id = vp_patients_mv.icustay_id; 

-- now make indicator variable if patient recieved fluids treatment.
select *, ( CASE
		WHEN cvvp_amount is not Null or mvvp_amount is not Null Then 1
    	ELSE 0
    	END ) AS vp_ind
into HE_ds_vp_trt
from he_ds_vp_all;

drop table he_ds_vpcv;
drop table he_ds_vp_all;

-- step 9. Make dataset that provides each patient's treatment or lack of treatment.
WITH he_fluid_temp AS (
						select icustay_id as icd, iv_fluid_ind
    					from he_ds_fluids_trt
)
 , he_trt_temp AS (
    select icustay_id, he_starttime, he_endtime, admittime, he_duration,
           pre_se_cr, post_se_cr, chg_ser, iv_fluid_ind, vp_ind
    from he_ds_vp_trt 
    inner join he_fluid_temp on he_ds_vp_trt.icustay_id = he_fluid_temp.icd
)
select * 
into he_ds_trt
from he_trt_temp;

--
select count(icustay_id)
from he_ds_trt
where iv_fluid_ind = 1 and vp_ind = 0;

select *, EXTRACT(epoch FROM he_duration)/60 as he_duration_mins 
into he_ds_trt_complete
from he_ds_trt;
