set search_path to mimiciii;

-- Important Items to examine.
select * from d_items
where itemid in (1125, 225312, 220052, 6926, 52, 6702, 456, 220181);

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
Values (60, 1, 60.01, 1);

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
    FROM chartevents_adult WHERE itemid in (456,52)
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
    FROM chartevents_adult WHERE itemid in (456,52) 
),

hypoexits2 as(
    select h.*,
    cast('EXIT' AS text) as status
    FROM hypoexits h WHERE valuenum is not null and valuenum between (select exit_hypo_threshold from dual) and 180
                                       and windowend between (select exit_hypo_threshold from dual) and 180 
    								   and icustay_id is not null
    								   and itemid = windowend_itemid and he_time != windowendct
    								   and (windowendct - he_time) <= '01:00:00'
),
--select * from hypoexits2 limit 500;

allevents as(
  select aa.* from (select * from hypoentries2 union select * from hypoexits2) aa order by subject_id, he_time
),
--select * from allevents limit 1000;

allevents2 as(
	select *,
    lag(icustay_id, 1) over (partition by subject_id, icustay_id order by he_time) as prev_icustay,
    lag(status, 1) over (partition by subject_id, icustay_id order by he_time) as prev_status,
    lag(itemid,1 ) over (partition by subject_id, icustay_id order by he_time) as prev_itemid
    from allevents
),
--select * from allevents2 limit 500;

allevents3 as(
        select distinct icustay_id from allevents2 where ( icustay_id = prev_icustay and status != prev_status and itemid = prev_itemid )
),
--select * from allevents3;

allevents31 as(
	select a.* from allevents2 a 
    inner join allevents3 b
    on a.icustay_id = b.icustay_id
),
--select * from allevents31 limit 500;

allevents35 as(
  select subject_id, icustay_id, itemid, valuenum, he_time, status, 
    lag(icustay_id, 1) over (partition by subject_id, icustay_id order by he_time) as prev_icustay,
    lag(status, 1) over (partition by subject_id, icustay_id order by he_time) as prev_status,
    lead(icustay_id, 1) over (partition by subject_id, icustay_id order by he_time) as next_icustay,
    lead(status, 1) over (partition by subject_id, icustay_id order by he_time) as next_status,
    LEAD(he_time, 1) over (partition by SUBJECT_ID, ICUSTAY_ID order by he_time) as NEXT_CHARTTIME
  from allevents31
),
--select count(*) from allevents35 limit 500; 

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
select h.*, (h.he_offset - h.he_onset) as he_length into he_set from hetable1 h where he_onset is not null order by h.subject_id, h.he_onset;

/*
This statement finds all hypotensive episodes in all icustays.
*/

select * from he_set where he_offset is not null;

with hypo_counts as(
	select subject_id, icustay_id, count(icustay_id) from he_set where he_offset is not null group by subject_id,icustay_id order by subject_id
)
, hypo_counts1 as(
select * from hypo_counts where count=1
),
he_subset as(
	select * from he_set where icustay_id in ( select icustay_id from hypo_counts1 ) 
   ),
--select * from he_subset;
/* Pick only first HE for each subject */
usable_he0 AS(
	SELECT distinct h.subject_id, h.icustay_id, h.itemid,
   	   first_value(h.he_onset) OVER (PARTITION BY h.subject_id, h.icustay_id ORDER BY h.he_onset ) as he_onset_first,
       first_value( h.he_offset ) OVER (PARTITION BY h.subject_id, h.icustay_id ORDER BY case when h.he_offset is not null then 0 else 1 end ASC, h.he_onset 
       rows unbounded preceding )  as he_offset_first 
     FROM he_subset h
),
usable_he_05 as(
	select subject_id, count(subject_id) as count_s, count(icustay_id) as count_i from usable_he0 group by subject_id having count(icustay_id) = 1 order by subject_id
)
select * from usable_he0 where subject_id in (select subject_id from usable_he_05 ) order by subject_id;

select distinct on (subject_id) * from usable_he0 where he_offset is not null
usable_he01 as (
	select * from usable_he0 where he_onset is not null and he_offset is not null
    ),
usable_he0_5 as(
	select distinct on (subject_id, he_onset) * from usable_he01 order by subject_id 
    )
select * --into he_set_final 
from usable_he0_5 where subject_id in (select subject_id from usable_he0_5 group by subject_id having count(subject_id)=1) 
order by subject_id;



/*
usable_he1 as(
    select *, (he_offset - he_onset) as he_length 
    from usable_he0_5 where he_onset is not null and he_offset is not null 
),
usable_he2 as(
select distinct * from usable_he1
)
*/

--select * from he_set_final;
--drop table he_set_final;

select k, percentile_disc(k) within group (order by he_length)
from he_set_final, generate_series(0.25, 0.75, 0.25) as k 
--where dbsource='carevue'
group by k;

select itemid, count(subject_id) from he_set_final group by itemid;
select * from he_set_final;
/*
This statement finds all icustays with one hypotensive episodes.
*/
with hh as(
select * --into he_set1 
from HE_SET
where subject_id in (
                        select subject_id from HE_SET_final
                        group by subject_id
                        having count(icustay_id) = 1 and count(subject_id)=1
)
 )
select itemid, count( distinct subject_id) from hh group by itemid;

-- intime is admittime to icu for each patient.
select dbsource, hadm_id, HE_set_final.*, intime, los
	into HE_Cohort
    from he_set_final
	INNER JOIN icustays
	ON HE_set_final.icustay_id = icustays.icustay_id 
	order by HE_set_final.icustay_id ;
    
/* Find quartiles for he_length to make sure it matches with "Interrogating a clinical database to study treatment of hypotension in the critically ill" */
select k, percentile_disc(k) within group (order by he_length)
from he_cohort, generate_series(0.25, 0.75, 0.25) as k 
--where dbsource='carevue'
group by k;

select * from he_cohort order by he_length;


--drop table HE_Cohort; 
--drop table HE_SET;
--drop table HE_SET1;
--drop table hypo_cohort_final_cv;
--drop table dual;

-- Remove pateints with CMO's.
WITH icustay_CMO AS (
    select h.icustay_id, chartevents_adult.itemid, chartevents_adult.value, chartevents_adult.charttime as CMO_time, h.he_length, h.he_onset, h.he_offset 
         from chartevents_adult
         inner join he_cohort h on chartevents_adult.icustay_id = h.icustay_id
          where  value in ('Comfort Measures','Comfort measures only','CMO') and
                        (h.he_onset - chartevents_adult.charttime ) between '-24:00:00' and '24:00:00' or
                   value in ('Comfort Measures', 'Comfort measures only','CMO') and
                        (h.he_offset - chartevents_adult.charttime ) between '-24:00:00' and '24:00:00'
    		)
select * into hypo_cohort_final_cv from he_cohort
where icustay_id not in (select distinct icustay_id from icustay_CMO);

/* Find quartiles for he_length to make sure it matches with "Interrogating a clinical database to study treatment of hypotension in the critically ill" */
select k, percentile_disc(k) within group (order by he_length)
from hypo_cohort_final_cv, generate_series(0.25, 0.75, 0.25) as k 
--where dbsource='carevue'
group by k;

select itemid, count(subject_id) from  he_set_final group by itemid;

--select * from hypo_cohort_final_cv;
--drop table hypo_cohort_final_cv;
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


