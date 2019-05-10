set search_path to mimiciii;

-- Important Items to examine.
select * from d_items
where itemid in (1125, 225312, 220052, 6926, 52, 6702, 456, 220181);

select * from chartevents_adult limit 5;
select * from inputevents_cv limit 5;
select * from inputevents_mv limit 5;


/*

/* Remove patients younger than 15. */
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
    subject_id, dob, gender, first_admittime, 
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
	,sub_age AS (
    SELECT *
    FROM age
     )
select sub_age.*, icustays.icustay_id
into subject_age
from sub_age
inner join icustays on sub_age.subject_id = icustays.subject_id;
*/

/*
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
*/


/* 

  Thank you to Dr. Rithi Koshari for sharing some SQL code for extracting patients with hypotensive episodes in MIMIC-II
  
  This query returns DBSOURCE, HADM_ID, SUBJECT_ID, ICUSTAY_ID, HE_ONSET, HE_OFFSET, HE_LENGTH, LOS for 
  subjects in MIMIC-III DB on criteria:
  
  Entry: Time when first of 2 measurements of less than 60 mmHg within 1 hour was made
  Exit: Time when first of 2 measurements of greater than 60 mmHg within 1 hour was made   
    
****
  In this query, an EVENT is a record in the table with the status ENTRY or EXIT 
  being assigned to the beginnings of windows that qualify as ENTRY or EXIT 
  criteria as explained above. This term will be used from here on out.
****
  
Hypoentries: gets charttimes and values of windows of 1 hour that qualify as entry criteria
  
Hypoentries2: returns charttime of last measurement within 1 hour window. This table
			  finds all patients that fit the definition of the start of a hypotensive episode.

Hypoexits: gets charttimes and values of windows of 1 hour that qualify as exit criteria

Allevents: Gathers all ENTRY and EXIT events into one table, as well as
  information about the previous record for filtering
   
HEtable1/final query: Assembles an HE event if the next record has the same
  ICUSTAY_ID, the current event is ENTRY and the next event is EXIT.
  
After this query, the SQL file "Hypo_fluid_pressor_mimiciii.sql" will use hypo_cohort_final to obtain patients who recieved treatment.

*/

CREATE TABLE dual        
	(
     enter_hypo_threshold int, 
     enter_windownum int,
     exit_hypo_threshold int, 
     exit_windownum int
);
INSERT INTO dual
Values (60, 1, 60, 1);

-- drop table dual;
-- drop table he_set4;
-- drop table he_set;
-- drop table he_cohort;
-- select * from dual;

/* BP measurements lower than 40 can be too low to be realistic sometimes. */

With hypoentries as (
select subject_id, icustay_id, itemid, valuenum, charttime as he_time,
    row_number()  over (partition by subject_id, icustay_id order by charttime) as rn,
    last_value(charttime) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowendct,
    last_value(valuenum) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowend,
    last_value(itemid) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowend_itemid
    FROM chartevents_adult WHERE itemid in (52, 456)
),
--select * from hypoentries limit 50;

hypoentries2 as (
select h.*,
 cast('ENTRY' AS text) as status
 FROM hypoentries h WHERE valuenum is not null and valuenum between 40 and (select enter_hypo_threshold from dual) 
                                      and windowend between 40 and (select enter_hypo_threshold from dual) and icustay_id is not null
    								  and itemid = windowend_itemid and he_time != windowendct
    								  and (windowendct - he_time) <= '01:00:00'
),
--select * from hypoentries2 limit 500;

/* BP measurements greater than 180 can be too low to be realistic sometimes. */

hypoexits as(
	select subject_id, icustay_id, itemid, valuenum, charttime as he_time,
    row_number()  over (partition by subject_id, icustay_id order by charttime) as rn,
    last_value(charttime) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowendct,
    last_value(valuenum) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowend,
    last_value(itemid) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowend_itemid
    FROM chartevents_adult WHERE itemid in (52, 456) 
),

hypoexits2 as(
    select h.*,
    cast('EXIT' AS text) as status
    FROM hypoexits h WHERE valuenum is not null and valuenum > (select exit_hypo_threshold from dual) and valuenum <= 180
                                       and windowend > (select exit_hypo_threshold from dual) and windowend <= 180 
    								   and icustay_id is not null
    								   and itemid = windowend_itemid and he_time != windowendct
    								   and (windowendct - he_time) <= '01:00:00'
),
--select * from hypoexits2 limit 500;

allevents as(
  select aa.* from (select * from hypoentries2 union select * from hypoexits2) aa order by subject_id, he_time
)
--select * from allevents limit 10000;

, allevents2 as(
	select *,
    lag(icustay_id, 1) over (partition by subject_id, icustay_id order by he_time) as prev_icustay,
    lag(status, 1) over (partition by subject_id, icustay_id order by he_time) as prev_status,
    lag(itemid,1 ) over (partition by subject_id, icustay_id order by he_time) as prev_itemid
    from allevents
)
--select * from allevents2 limit 10000;
, allevents3_icu as(
        select icustay_id from allevents2 where ( icustay_id = prev_icustay and status != prev_status and itemid = prev_itemid ) 
    								
)
, hypo_events as(
	select * from allevents where icustay_id in (select * from allevents3_icu)
 )
--select * from hypo_events limit 5000;    

, allevents35 as(
  select subject_id, icustay_id, itemid, valuenum, he_time, status, 
    lag(icustay_id, 1) over (partition by subject_id, icustay_id order by he_time) as prev_icustay,
    lag(status, 1) over (partition by subject_id, icustay_id order by he_time) as prev_status,
    lead(icustay_id, 1) over (partition by subject_id, icustay_id order by he_time) as next_icustay,
    lead(status, 1) over (partition by subject_id, icustay_id order by he_time) as next_status,
    LEAD(he_time, 1) over (partition by SUBJECT_ID, ICUSTAY_ID order by he_time) as NEXT_CT
  from hypo_events
)
--select * from allevents35 limit 10000; 

, hetable1 as(
  select SUBJECT_ID, ICUSTAY_ID, itemid,
    case 
      when (status = 'ENTRY') 
        then he_time else null 
    end he_onset,
    case 
      when (status = 'ENTRY' and ICUSTAY_ID = NEXT_ICUstay and NEXT_STATUS = 'EXIT')
        then NEXT_CT else null 
    end he_offset
  from allevents35
)
select h.*, (h.he_offset - h.he_onset) as he_length into he_set 
from hetable1 h where he_onset is not null order by h.subject_id, h.he_onset;

/*
The following statements find subjects with one HE.
*/

with numHE as (
	select * from he_set where he_offset is not null
)
, count_subject as (
	select subject_id, count(subject_id) as count_id from numHE group by subject_id having count(subject_id)=1
)
, he_set1 as(
	select * from he_set where subject_id in (select subject_id from count_subject)
)
, usable_he AS(
SELECT distinct h.subject_id, h.icustay_id, h.itemid,
   	   first_value(h.he_onset) OVER (PARTITION BY h.subject_id, h.icustay_id ORDER BY h.he_onset ) he_onset,
       first_value( h.he_offset ) OVER (PARTITION BY h.subject_id, h.icustay_id ORDER BY case when h.he_offset is not null then 0 else 1 end ASC, h.he_offset 
     rows unbounded preceding ) he_offset
    FROM he_set1 h
)
, usable_he2 AS(
    select *, he_offset - he_onset as he_length from usable_he order by subject_id
 )
 -- Some subjects have BP recorded under 52 and 456 at the same times and these duplicates need to be reduced to one HE.
, distinct_HE AS(
	select distinct on (subject_id) * from usable_he2 where subject_id in (select subject_id from usable_he2 group by subject_id having count(subject_id)>1) and he_offset is not null
)
, hypo_cohort_union AS(
    select * from usable_he2 where subject_id not in (select subject_id from distinct_HE) and he_offset is not null
    union
	select * from distinct_HE
)
select * into hypo_cohort from hypo_cohort_union order by subject_id;

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

with usable_he AS(
SELECT distinct h.subject_id, icustay_id, first_value(h.icustay_id) OVER (PARTITION BY h.subject_id ORDER BY h.he_onset) icustay_id_first,
       first_value(h.itemid) OVER (PARTITION BY h.subject_id ORDER BY h.he_onset) first_itemid
--   	   first_value(h.he_onset) OVER (PARTITION BY h.subject_id, h.icustay_id ORDER BY h.he_onset ) he_onset,
--       first_value( h.he_offset ) OVER (PARTITION BY h.subject_id, h.icustay_id ORDER BY case when h.he_offset is not null then 0 else 1 end ASC, h.he_offset 
--       rows unbounded preceding ) he_offset
    FROM he_set h where he_offset is not null
)
--select * from usable_he;
, first_icu AS(
	select h.* from he_set h, usable_he u where h.icustay_id = u.icustay_id_first and h.itemid = u.first_itemid order by h.subject_id
)
--select * from first_icu;
, he_eps AS(
	select *, he_offset - he_onset as he_length from first_icu where he_offset - he_onset is not null
)
, he_eps1 AS(
	select subject_id, count(subject_id) from he_eps group by subject_id having count(subject_id) = 1 order by subject_id
)
--select * from he_eps1;
, usable_he2 AS(
SELECT distinct h.subject_id, icustay_id, 
         first_value(h.itemid) OVER (PARTITION BY h.subject_id ORDER BY h.he_onset) itemid,
    	 first_value(h.he_onset) OVER (PARTITION BY h.subject_id, h.icustay_id ORDER BY h.he_onset ) he_onset,
         first_value( h.he_offset ) OVER (PARTITION BY h.subject_id, h.icustay_id ORDER BY case when h.he_offset is not null then 0 else 1 end ASC, h.he_offset 
         rows unbounded preceding ) he_offset
    FROM first_icu h where he_offset is not null
)
select *, he_offset - he_onset as he_length into hypo_cohort_icu 
from usable_he2 where subject_id in (select subject_id from he_eps1) order by subject_id;


drop table hypo_cohort_icu;

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

with usable_he AS(
SELECT distinct h.subject_id, icustay_id, 
       dense_rank() over (partition by h.subject_id order by h.icustay_id)  + dense_rank() over (partition by h.subject_id order by h.icustay_id desc) -1 as count_icu,
       dense_rank() over (partition by h.subject_id, h.icustay_id order by h.he_onset)  +  dense_rank() over (partition by h.subject_id, h.icustay_id order by h.he_onset desc) -1 
       as count_events,
       dense_rank() over (partition by h.subject_id, h.icustay_id order by h.itemid)  +  dense_rank() over (partition by h.subject_id, h.icustay_id order by h.itemid desc) -1 
       as count_itemid,
       first_value(h.he_offset) OVER (PARTITION BY h.subject_id ORDER BY h.he_offset) first_offset,
       first_value(h.itemid) OVER (PARTITION BY h.subject_id ORDER BY h.he_onset) first_item
       FROM he_set h where he_offset is not null
)
--select * from usable_he order by subject_id;
, he_subjs AS(
    select h.* from he_set h, usable_he u where (count_events = 1 and count_icu = 1) 
    and h.he_offset=u.first_offset 
    and he_offset is not null 
    order by subject_id
)
, usable_he2 AS(
SELECT distinct h.subject_id, icustay_id, itemid,
    	 first_value(h.he_onset) OVER (PARTITION BY h.subject_id, h.icustay_id ORDER BY h.he_onset ) he_onset,
         first_value( h.he_offset ) OVER (PARTITION BY h.subject_id, h.icustay_id ORDER BY case when h.he_offset is not null then 0 else 1 end ASC, h.he_offset 
         rows unbounded preceding ) he_offset
    FROM he_subjs h where he_offset is not null
)
select *, he_offset - he_onset as he_length into hypo_cohort_icu 
from usable_he2
order by subject_id;

/* This is the code to use above */
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
, one_stay AS(
	select * from usable_he where count_icu = 1 order by subject_id
)
--select * from one_stay;
, first_icu AS(
	select h.* from he_set h, one_stay o where h.icustay_id = o.icustay_id and h.itemid = o.first_itemid and he_offset is not null order by h.subject_id
)
--select * from first_icu;
, he_eps AS(
	select *, he_offset - he_onset as he_length from first_icu where he_offset - he_onset is not null
)
--select * from he_eps;
, he_eps1 AS(
	select subject_id, count(subject_id) from he_eps group by subject_id having count(subject_id) = 1 order by subject_id
)
, usable_he2 AS(
SELECT distinct h.subject_id, icustay_id, itemid,
         --first_value(h.itemid) OVER (PARTITION BY h.subject_id ORDER BY h.he_onset) itemid,
    	 first_value(h.he_onset) OVER (PARTITION BY h.subject_id, h.icustay_id ORDER BY h.he_onset ) he_onset,
         first_value( h.he_offset ) OVER (PARTITION BY h.subject_id, h.icustay_id ORDER BY case when h.he_offset is not null then 0 else 1 end ASC, h.he_offset 
         rows unbounded preceding ) he_offset
    FROM first_icu h where he_offset is not null
)
select *, he_offset - he_onset as he_length --into hypo_cohort_icu 
from usable_he2 where subject_id in (select subject_id from he_eps1) order by subject_id;


drop table hypo_cohort_icu;

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/* New approach with each subject with one ICU Stay */


-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
WITH usable_subj AS(
	SELECT h.subject_id
       FROM he_set h where (he_offset - he_onset) is not null group by h.subject_id having count( subject_id ) = 1
)
--select * from usable_subj order by subject_id; -- 3128
, usable_he AS(
	SELECT distinct h.subject_id, icustay_id, itemid, he_onset, he_offset,
         first_value(h.itemid) OVER (PARTITION BY h.subject_id ORDER BY h.he_onset) itemid_1,
         first_value(h.icustay_id) OVER (PARTITION BY h.subject_id ORDER BY h.he_onset) icustay_1,
    	 first_value(h.he_onset) OVER (PARTITION BY h.subject_id, h.icustay_id ORDER BY h.he_onset ) he_onset_1,
         first_value( h.he_offset ) OVER (PARTITION BY h.subject_id, h.icustay_id ORDER BY case when h.he_offset is not null then 0 else 1 end ASC, h.he_offset 
         rows unbounded preceding ) he_offset_1
    FROM he_set h where he_offset is not null
)
select distinct subject_id from usable_he 
where he_onset = he_onset_1 and subject_id in (select subject_id from usable_he where he_onset = he_onset_1 and icustay_id = icustay_1 )
order by subject_id;

select subject_id from usable_he group by subject_id having count(subject_id)=4 and count(distinct itemid)=2 order by subject_id;



, one_stay AS(
	select * from usable_he where count_icu = 1 order by subject_id
)
--select * from one_stay;
, first_icu AS(
	select h.* from he_set h, one_stay o where h.icustay_id = o.icustay_id and h.itemid = o.first_itemid order by h.subject_id
)
--select * from first_icu;
, he_eps AS(
	select *, he_offset - he_onset as he_length from first_icu where he_offset - he_onset is not null
)
--select * from he_eps;
, he_eps1 AS(
	select subject_id, count(subject_id) from he_eps group by subject_id having count(subject_id) = 1 order by subject_id
)
, usable_he2 AS(
SELECT distinct h.subject_id, icustay_id, itemid,
         first_value(h.itemid) OVER (PARTITION BY h.subject_id ORDER BY h.he_onset) itemid,
    	 first_value(h.he_onset) OVER (PARTITION BY h.subject_id, h.icustay_id ORDER BY h.he_onset ) he_onset,
         first_value( h.he_offset ) OVER (PARTITION BY h.subject_id, h.icustay_id ORDER BY case when h.he_offset is not null then 0 else 1 end ASC, h.he_offset 
         rows unbounded preceding ) he_offset
    FROM first_icu h where he_offset is not null
)
select *, he_offset - he_onset as he_length --into hypo_cohort_icu 
from usable_he2 where subject_id in (select subject_id from he_eps1) order by subject_id;


drop table hypo_cohort_icu;

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------


/*
This statement finds the number of hypotensive episodes per icustay and finds patients with unique HE's at first stay.
*/

-- Find all subjects with multiple ICUStays but one HE in first ICUStay.
with numHE as (
	select * from he_set where he_offset is not null
)
, count_icustay as (
	select subject_id, count(distinct icustay_id) as count_id from numHE group by subject_id having count(distinct icustay_id)>1
)
, he_set1 as(
	select * from he_set where subject_id in (select subject_id from count_icustay)
)
/* Pick only first HE for each subject */
, usable_he0 AS(
	SELECT *, first_value(h.icustay_id) OVER (PARTITION BY h.subject_id ORDER BY h.he_onset) icustay_id_first
	FROM he_set1 h
)
, usable_he1 AS(
	select subject_id, icustay_id, itemid, he_onset, he_offset, he_length from usable_he0 where icustay_id = icustay_id_first order by subject_id
)
--select * from usable_he1;
/* Find Length of Episode by Finding the Start and End of HE */
, usable_he2 AS(
SELECT distinct h.subject_id, h.icustay_id, h.itemid,
   	   first_value(h.he_onset) OVER (PARTITION BY h.subject_id, h.icustay_id ORDER BY h.he_onset ) he_onset,
       first_value( h.he_offset ) OVER (PARTITION BY h.subject_id, h.icustay_id ORDER BY case when h.he_offset is not null then 0 else 1 end ASC, h.he_offset 
     rows unbounded preceding ) he_offset
    FROM usable_he1 h
)
, usable_he25 AS(
    select subject_id from usable_he2 group by subject_id having count(subject_id) = 2
)
--select * from usable_he2 order by subject_id;
, usable_he3 AS(
	select *, (h.he_offset - h.he_onset) as he_length from usable_he2 h where subject_id not in (select subject_id from usable_he25 ) order by subject_id
)
, all_he as(
    select * from usable_he3
    --union 
    --select * from hypo_cohort order by subject_id
)
select * from all_he;

/*
/* Pick only first HE for each subject with multiple ICU Stays */

with usable_he0 AS(
	SELECT *, first_value(h.icustay_id) OVER (PARTITION BY h.subject_id ORDER BY h.he_onset) icustay_id_first
    FROM he_set1 h
)
, usable_he1 AS(
	select * from usable_he0 where icustay_id = icustay_id_first order by subject_id
)
, hypo_counts_subj1 as(
	select subject_id, count(subject_id) as count_subj from usable_he1 group by subject_id order by subject_id
)
select * into hypo_cohort from usable_he1 where subject_id in (select subject_id from hypo_counts_subj1 where count_subj=1);

select k, percentile_disc(k) within group (order by he_length)
from hypo_cohort, generate_series(0.25, 0.75, 0.25) as k 
--where dbsource='carevue'
group by k;
*/

/* If you want to drop any tables

--drop table hypo_cohort; 
--drop table HE_SET;
--drop table HE_SET1;
--drop table hypo_cohort_final_cv;
--drop table dual;

*/

-- Remove pateints with CMO's.
WITH icustay_CMO AS (
    select h.icustay_id, chartevents_adult.itemid, chartevents_adult.value, chartevents_adult.charttime as CMO_time, h.he_length, h.he_onset, h.he_offset 
         from chartevents_adult
         inner join hypo_cohort_icu h on chartevents_adult.icustay_id = h.icustay_id
          where  value in ('Comfort Measures','Comfort measures only','CMO') and
                        (h.he_onset - chartevents_adult.charttime ) between '-24:00:00' and '24:00:00' or
                   value in ('Comfort Measures', 'Comfort measures only','CMO') and
                        (h.he_offset - chartevents_adult.charttime ) between '-24:00:00' and '24:00:00'
    		)
select * into hypo_cohort_final_cv 
from hypo_cohort_icu
where icustay_id not in (select icustay_id from icustay_CMO) order by he_length;

-- intime is admittime to icu for each patient.
with icu_hypo AS(
    select dbsource, hadm_id, h.*, intime, los
        from hypo_cohort_final_cv  h
        INNER JOIN icustays
        ON h.icustay_id = icustays.icustay_id 
        order by h.icustay_id
)
, pats_hypo AS(
 select h.*, gender, dob 
 from icu_hypo h
 INNER JOIN patients p
 ON p.subject_id = h.subject_id
 order by h.subject_id
)
, inhm_hypo AS(
  select p.*, HOSPITAL_EXPIRE_FLAG
  from pats_hypo p
  inner join admissions a
  ON p.hadm_id = a.hadm_id
)
select *, age(he_offset , dob) as age into hypo_cohort_cv 
from inhm_hypo;

drop table hypo_cohort_cv; 

select * from hypo_cohort_cv where age >= '18 years';
select * from hypo_cohort_cv where gender='M';
select count(*) from hypo_cohort_cv where HOSPITAL_EXPIRE_FLAG=1; --507
select count(*) from hypo_cohort_cv; --3576


/* Find quartiles for he_length to make sure it matches with "Interrogating a clinical database to study treatment of hypotension in the critically ill" */
select k, percentile_disc(k) within group (order by los)
from hypo_cohort_cv, generate_series(0.25, 0.75, 0.25) as k 
group by k;

select itemid, count(subject_id) from hypo_cohort_cv group by itemid;

--select * from hypo_cohort_final_cv;
--drop table hypo_cohort_final_cv;
--drop table hypo_cohort_cv;
--select avg(los) from hypo_cohort_final;

/*

Extra Tables that can be removed:

-- drop table he_set1_cmo;
-- drop table he_set;
-- drop table he_set1;
-- drop table he_cohort;

*/

---------------------------------------------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------------------------------------------

drop table he_set;
drop table he_set1;
drop table he_cohort;

/*

Next, find patients suffering from hypotensive episodes in the metavision database. 

*/

With hypoentries as (
select subject_id, icustay_id, itemid, valuenum, charttime as he_time,
    row_number()  over (partition by subject_id, icustay_id order by charttime) as rn,
    last_value(charttime) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowendct,
    last_value(valuenum) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowend,
    last_value(itemid) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowend_itemid
    FROM chartevents_adult WHERE itemid = 52 --in ( 220181, 220052, 225312 )
),
--select * from hypoentries limit 50;

hypoentries2 as (
select h.*,
 cast('ENTRY' AS text) as status
 FROM hypoentries h WHERE valuenum is not null and valuenum between 40 and (select enter_hypo_threshold from dual) 
                                      and windowend between 40 and (select enter_hypo_threshold from dual)and icustay_id is not null
    								  and itemid = windowend_itemid and he_time != windowendct
),
--select * from hypoentries2 limit 500;

hypoexits as(
select subject_id, icustay_id, itemid, valuenum, charttime as he_time,
    row_number()  over (partition by subject_id, icustay_id order by charttime) as rn,
    last_value(charttime) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowendct,
    last_value(valuenum) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowend,
    last_value(itemid) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowend_itemid
    FROM chartevents_adult WHERE itemid = 52 --in ( 220181, 220052, 225312 )
),

hypoexits2 as(
    select h.*,
    cast('EXIT' AS text) as status
    FROM hypoexits h WHERE valuenum is not null and valuenum > (select exit_hypo_threshold from dual) and valuenum <= 180
                                       and windowend > (select exit_hypo_threshold from dual) and valuenum <= 180 
    								   and icustay_id is not null
    								   and itemid = windowend_itemid and he_time != windowendct
),
--select * from hypoexits2 limit 500;

allevents as(
  select * from (select * from hypoentries2 union select * from hypoexits2) aa order by subject_id, he_time
),
--select * from allevents limit 1000;

allevents2 as(
	select *,
    lag(icustay_id, 1) over (partition by subject_id, icustay_id order by he_time) as prev_icustay,
    lag(status, 1) over (partition by subject_id, icustay_id order by he_time) as prev_status,
    lag(itemid, 1) over (partition by subject_id, icustay_id order by he_time) as prev_itemid
    from allevents
),
--select * from allevents2 limit 500;

allevents3 as(
        select distinct icustay_id from allevents2 where ( icustay_id = prev_icustay and status != prev_status and itemid = prev_itemid )
),
--select * from allevents3;
--6,093 subjects

allevents31 as(
	select a.* from allevents2 a 
    inner join allevents3 b
    on a.icustay_id = b.icustay_id
),
--select * from allevents31;

allevents35 as(
  select subject_id, icustay_id, itemid, valuenum, he_time, status,
    lag(icustay_id, 1) over (partition by subject_id, icustay_id order by he_time) as prev_icustay,
    lag(status, 1) over (partition by subject_id, icustay_id order by he_time) as prev_status,
    lead(icustay_id, 1) over (partition by subject_id, icustay_id order by he_time) as next_icustay,
    lead(status, 1) over (partition by subject_id, icustay_id order by he_time) as next_status,
    LEAD(he_time, 1) over (partition by SUBJECT_ID, ICUSTAY_ID order by he_time) as NEXT_CHARTTIME
  from allevents31
),
--select * from allevents35;

hetable1 as(
  select SUBJECT_ID, ICUSTAY_ID, itemid,
    case 
      when (status = 'ENTRY') 
        then he_time else null 
    end he_onset,
    case 
      when (status = 'ENTRY' and ICUSTAY_ID = NEXT_ICUstay and NEXT_STATUS = 'EXIT')
        then NEXT_CHARTTIME else null 
    end he_offset
  from allevents35
)
--select h.* from hetable1 h;
select h.*, (h.he_offset - h.he_onset) as he_length into he_set_52 from hetable1 h where he_onset is not null or he_offset is not null order by h.subject_id, h.he_onset;

/* Pick only first HE for each subject */
with usable_he0 AS(
SELECT h.subject_id, h.icustay_id, h.itemid,
   	   first_value(h.he_onset) OVER (PARTITION BY h.subject_id, h.icustay_id ORDER BY h.he_onset ) he_onset,
       first_value( h.he_offset ) OVER (PARTITION BY h.subject_id, h.icustay_id ORDER BY case when h.he_offset is not null then 0 else 1 end ASC, h.he_onset 
     rows unbounded preceding ) he_offset FROM he_set_52 h
),
--select distinct * from usable_he0 order by subject_id;
usable_he1 as(
    select *, (he_offset - he_onset) as he_length 
    from usable_he0 where he_onset is not null and he_offset is not null 
),
usable_he2 as(
select distinct * from usable_he1
)
select * into he_set_52_final 
from usable_he2 where subject_id in (select subject_id from usable_he2 group by subject_id having count(subject_id)=1) 
order by subject_id;

/*
This statement finds all hypotensive episodes in all icustays.
*/

--select h.*, (h.he_offset - h.he_onset) as he_length into HE_SET2 from hetable1 h where he_onset is not null and he_offset is not null order by h.subject_id, h.he_onset;

/*
This statement finds all icustays with one hypotensive episodes.
*/

select * into he_set3 from HE_SET2 
where subject_id in (
                        select subject_id from HE_SET2
                        group by subject_id
                        having count(icustay_id) = 1
);

select dbsource, hadm_id, HE_set3.*, intime as admittime, los
	into HE_Cohort2
    from he_set3
	INNER JOIN icustays
	ON HE_set3.icustay_id = icustays.icustay_id 
	order by HE_set3.icustay_id ;
 
select k, percentile_disc(k) within group (order by los)
from he_cohort2, generate_series(0.25, 0.75, 0.25) as k 
--where dbsource='carevue'
group by k;

select subject_id, count(subject_id) from he_set2 group by subject_id having count(subject_id)=1 order by subject_id;
select count( distinct icustay_id) from he_set2;
select count( distinct icustay_id) from he_set1;


--drop table HE_SET2;
--drop table HE_SET3;
--drop table HE_Cohort2; 
--drop table dual;
-- drop table hypo_cohort_final_mv;

--------------------------------------------------------
With hypoentries as (
select subject_id, icustay_id, itemid, valuenum, charttime as he_time,
    row_number()  over (partition by subject_id, icustay_id order by charttime) as rn,
    last_value(charttime) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowendct,
    last_value(valuenum) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowend,
    last_value(itemid) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowend_itemid
    FROM chartevents_adult WHERE itemid = 6702 --in ( 220181, 220052, 225312 )
),
--select * from hypoentries limit 50;

hypoentries2 as (
select h.*,
 cast('ENTRY' AS text) as status
 FROM hypoentries h WHERE valuenum is not null and valuenum between 40 and (select enter_hypo_threshold from dual) 
                                      and windowend between 40 and (select enter_hypo_threshold from dual)and icustay_id is not null
    								  and itemid = windowend_itemid and he_time != windowendct
),
--select * from hypoentries2 limit 500;

hypoexits as(
select subject_id, icustay_id, itemid, valuenum, charttime as he_time,
    row_number()  over (partition by subject_id, icustay_id order by charttime) as rn,
    last_value(charttime) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowendct,
    last_value(valuenum) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowend,
    last_value(itemid) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowend_itemid
    FROM chartevents_adult WHERE itemid = 6702 --in ( 220181, 220052, 225312 )
),

hypoexits2 as(
    select h.*,
    cast('EXIT' AS text) as status
    FROM hypoexits h WHERE valuenum is not null and valuenum > (select exit_hypo_threshold from dual) and valuenum <= 180
                                       and windowend > (select exit_hypo_threshold from dual) and valuenum <= 180 
    								   and icustay_id is not null
    								   and itemid = windowend_itemid and he_time != windowendct
),
--select * from hypoexits2 limit 500;

allevents as(
  select * from (select * from hypoentries2 union select * from hypoexits2) aa order by subject_id, he_time
),
--select * from allevents limit 1000;

allevents2 as(
	select *,
    lag(icustay_id, 1) over (partition by subject_id, icustay_id order by he_time) as prev_icustay,
    lag(status, 1) over (partition by subject_id, icustay_id order by he_time) as prev_status,
    lag(itemid, 1) over (partition by subject_id, icustay_id order by he_time) as prev_itemid
    from allevents
),
--select * from allevents2 limit 500;

allevents3 as(
        select distinct icustay_id from allevents2 where ( icustay_id = prev_icustay and status != prev_status and itemid = prev_itemid )
),
--select * from allevents3;
--6,093 subjects

allevents31 as(
	select a.* from allevents2 a 
    inner join allevents3 b
    on a.icustay_id = b.icustay_id
),
--select * from allevents31;

allevents35 as(
  select subject_id, icustay_id, itemid, valuenum, he_time, status,
    lag(icustay_id, 1) over (partition by subject_id, icustay_id order by he_time) as prev_icustay,
    lag(status, 1) over (partition by subject_id, icustay_id order by he_time) as prev_status,
    lead(icustay_id, 1) over (partition by subject_id, icustay_id order by he_time) as next_icustay,
    lead(status, 1) over (partition by subject_id, icustay_id order by he_time) as next_status,
    LEAD(he_time, 1) over (partition by SUBJECT_ID, ICUSTAY_ID order by he_time) as NEXT_CHARTTIME
  from allevents31
),
/*
allevents36 as(
    select * from allevents35 
    where status != prev_status or status != next_status
), 
*/
hetable1 as(
  select SUBJECT_ID, ICUSTAY_ID, itemid,
    case 
      when (status = 'ENTRY') 
        then he_time else null 
    end he_onset,
    case 
      when (status = 'ENTRY' and ICUSTAY_ID = NEXT_ICUstay and NEXT_STATUS = 'EXIT')
        then NEXT_CHARTTIME else null 
    end he_offset
  from allevents35
)
--select h.*, (h.he_offset - h.he_onset) as he_length from hetable1 h where he_onset is not null and he_offset is not null order by h.subject_id, h.he_onset;
--select h.*, (h.he_offset - h.he_onset) as he_length from hetable1 h where he_onset is not null and he_offset is not null order by h.subject_id, h.he_onset;
select h.*, (h.he_offset - h.he_onset) as he_length into he_set_6702 from hetable1 h where he_onset is not null or he_offset is not null order by h.subject_id, h.he_onset;

/* Pick only first HE for each subject */
with usable_he0 AS(
SELECT distinct h.subject_id, first_value(h.icustay_id) OVER (PARTITION BY h.subject_id ORDER BY h.he_onset) icustay_id, itemid,
  first_value(h.he_onset) OVER (PARTITION BY h.subject_id ORDER BY h.he_onset) he_onset,
  first_value(h.he_offset) OVER (PARTITION BY h.subject_id ORDER BY h.he_onset) he_offset FROM he_set_6702 h
)
select *, (he_offset - he_onset) as he_length into he_set_6702_final from usable_he0 where he_onset is not null and he_offset is not null;


--select h.*, (h.he_offset - h.he_onset) as he_length into HE_SET3 from hetable1 h where he_onset is not null and he_offset is not null order by h.subject_id, h.he_onset;

-------------------------------------
With hypoentries as (
select subject_id, icustay_id, itemid, valuenum, charttime as he_time,
    row_number()  over (partition by subject_id, icustay_id order by charttime) as rn,
    last_value(charttime) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowendct,
    last_value(valuenum) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowend,
    last_value(itemid) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowend_itemid
    FROM chartevents_adult WHERE itemid = 443 --in ( 220181, 220052, 225312 )
),
--select * from hypoentries limit 50;

hypoentries2 as (
select h.*,
 cast('ENTRY' AS text) as status
 FROM hypoentries h WHERE valuenum is not null and valuenum between 40 and (select enter_hypo_threshold from dual) 
                                      and windowend between 40 and (select enter_hypo_threshold from dual)and icustay_id is not null
    								  and itemid = windowend_itemid and he_time != windowendct
),
--select * from hypoentries2 limit 500;

hypoexits as(
select subject_id, icustay_id, itemid, valuenum, charttime as he_time,
    row_number()  over (partition by subject_id, icustay_id order by charttime) as rn,
    last_value(charttime) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowendct,
    last_value(valuenum) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowend,
    last_value(itemid) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowend_itemid
    FROM chartevents_adult WHERE itemid = 443 --in ( 220181, 220052, 225312 )
),

hypoexits2 as(
    select h.*,
    cast('EXIT' AS text) as status
    FROM hypoexits h WHERE valuenum is not null and valuenum > (select exit_hypo_threshold from dual) and valuenum <= 180
                                       and windowend > (select exit_hypo_threshold from dual) and valuenum <= 180 
    								   and icustay_id is not null
    								   and itemid = windowend_itemid and he_time != windowendct
),
--select * from hypoexits2 limit 500;

allevents as(
  select * from (select * from hypoentries2 union select * from hypoexits2) aa order by subject_id, he_time
),
--select * from allevents limit 1000;

allevents2 as(
	select *,
    lag(icustay_id, 1) over (partition by subject_id, icustay_id order by he_time) as prev_icustay,
    lag(status, 1) over (partition by subject_id, icustay_id order by he_time) as prev_status,
    lag(itemid, 1) over (partition by subject_id, icustay_id order by he_time) as prev_itemid
    from allevents
),
--select * from allevents2 limit 500;

allevents3 as(
        select distinct icustay_id from allevents2 where ( status != prev_status and itemid = prev_itemid )
),
--select * from allevents3;
--6,093 subjects

allevents31 as(
	select a.* from allevents2 a 
    inner join allevents3 b
    on a.icustay_id = b.icustay_id
),
--select * from allevents31;

allevents35 as(
  select subject_id, icustay_id, itemid, valuenum, he_time, status,
    lag(icustay_id, 1) over (partition by subject_id, icustay_id order by he_time) as prev_icustay,
    lag(status, 1) over (partition by subject_id, icustay_id order by he_time) as prev_status,
    lead(icustay_id, 1) over (partition by subject_id, icustay_id order by he_time) as next_icustay,
    lead(status, 1) over (partition by subject_id, icustay_id order by he_time) as next_status,
    LEAD(he_time, 1) over (partition by SUBJECT_ID, ICUSTAY_ID order by he_time) as NEXT_CHARTTIME
  from allevents31
),
--select * from allevents35;
/*
allevents36 as(
    select * from allevents35 
    where status != prev_status or status != next_status
), 
*/

hetable1 as(
  select SUBJECT_ID, ICUSTAY_ID, itemid,
    case 
      when (status = 'ENTRY') 
        then he_time else null 
    end he_onset,
    case 
      when (status = 'ENTRY' and ICUSTAY_ID = NEXT_ICUstay and NEXT_STATUS = 'EXIT')
        then NEXT_CHARTTIME else null 
    end he_offset
  from allevents35
)

--select h.*, (h.he_offset - h.he_onset) as he_length into HE_SET4 from hetable1 h where he_onset is not null and he_offset is not null order by h.subject_id, h.he_onset;
select h.*, (h.he_offset - h.he_onset) as he_length into he_set_443 from hetable1 h where he_onset is not null or he_offset is not null order by h.subject_id, h.he_onset;

/* Pick only first HE for each subject */
with usable_he0 AS(
SELECT distinct h.subject_id, first_value(h.icustay_id) OVER (PARTITION BY h.subject_id ORDER BY h.he_onset) icustay_id, itemid,
  first_value(h.he_onset) OVER (PARTITION BY h.subject_id ORDER BY h.he_onset) he_onset,
  first_value(h.he_offset) OVER (PARTITION BY h.subject_id ORDER BY h.he_onset) he_offset FROM he_set_443 h
)
select *, (he_offset - he_onset) as he_length into he_set_443_final from usable_he0 where he_onset is not null and he_offset is not null;

------
with hh as(
select * from he_set 
union 
select * from he_set2
),
dh as(
select distinct subject_id from hh group by subject_id having count(subject_id) = 1
), 
hh2 as(
    select * from he_set3
    union 
    select * from he_set4
),
dh2 as(
	select distinct subject_id from hh2 group by subject_id having count(subject_id) = 1
),
allh as(
	select * from dh union select * from dh2
)
select distinct subject_id from allh;

select count(*) from he_set4; -- 23407 + 20998 + 225 + 18

-- drop table he_set_456_final;
-- drop table he_set_52_final;
-- drop table he_set_6702_final;
-- drop table he_set_443_final;

select h.* from he_set_52_final h
inner join he_set_456_final h2
on h.subject_id = h2.subject_id;

with he_set_final as (
    select * from he_set_456_final 
    union
    select * from he_set_52_final
--    union
--    select * from he_set_6702_final
--    union
--    select * from he_set_443_final
)
select subject_id from he_set_final group by subject_id having count(subject_id)=1 order by subject_id;
select * from he_set_final order by subject_id; --44648

he_set_count1 as (
	select * from he_set_final where subject_id in (select subject_id from he_set_final group by subject_id having count(subject_id)=1) order by subject_id
)
--select * from he_set_count1;

select dbsource, hadm_id, he_set_count1.*, intime as admittime, los
	into HE_Cohort4
    from he_set_count1
	INNER JOIN icustays
	ON he_set_count1.icustay_id = icustays.icustay_id 
	order by he_set_count1.icustay_id ;

--drop table he_cohort4;
-- Remove pateints with CMO's.
WITH icustay_CMO AS (
    select h.icustay_id, chartevents_adult.itemid, chartevents_adult.value, chartevents_adult.charttime as CMO_time, h.he_length, h.he_onset, h.he_offset 
         from chartevents_adult
         inner join he_cohort4 h on chartevents_adult.icustay_id = h.icustay_id
          where  value in ('Comfort Measures','Comfort measures only','CMO') and
                        (h.he_onset - chartevents_adult.charttime ) between '-24:00:00' and '24:00:00' or
                   value in ('Comfort Measures','Comfort measures only','CMO') and
                        (h.he_offset - chartevents_adult.charttime ) between '-24:00:00' and '24:00:00'
    		)
select * into hypo_cohort_final_cv_all from he_cohort4
where icustay_id not in (select distinct icustay_id from icustay_CMO);
 
select k, percentile_disc(k) within group (order by he_length)
from he_cohort4, generate_series(0.25, 0.75, 0.25) as k 
--where dbsource='carevue'
group by k;

--drop table he_set1;
--drop table he_set2;
--drop table he_set3;
--drop table he_set4;
--drop table hypo_cohort_final_cv_all;

/* Final cohort combines both these tables */
select count(*) from  hypo_cohort_final_mv --order by subject_id;

with he_final as (
    select * from hypo_cohort_final_cv 
    union
    select * from hypo_cohort_final_mv 
),
--select subject_id, count(subject_id) from he_final group by subject_id having count(subject_id)=2;
hh as(
	select * from he_final where subject_id in ( select subject_id from he_final group by subject_id having count(subject_id)=1)
),
--select * from hh;

hhh as (
	select * from hh 
    )
select k, percentile_disc(k) within group (order by los)
from hhh, generate_series(0.25, 0.75, 0.25) as k 
--where dbsource='carevue'
group by k;




select * into hypo_cohort_final from he_final;

select * from hypo_cohort_final_mv;


