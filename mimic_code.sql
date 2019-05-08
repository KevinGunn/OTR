
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
    FROM hypoexits h WHERE valuenum is not null and valuenum between (select exit_hypo_threshold from dual) and 180
                                       and windowend between (select exit_hypo_threshold from dual) and 180 
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

, allevents3 as(
        select * from allevents2 where ( icustay_id = prev_icustay and status != prev_status and itemid = prev_itemid )
)
--select * from allevents3 limit 10000;

, allevents35 as(
  select subject_id, icustay_id, itemid, valuenum, he_time, status, 
    lag(icustay_id, 1) over (partition by subject_id, icustay_id order by he_time) as prev_icustay,
    lag(status, 1) over (partition by subject_id, icustay_id order by he_time) as prev_status,
    lead(icustay_id, 1) over (partition by subject_id, icustay_id order by he_time) as next_icustay,
    lead(status, 1) over (partition by subject_id, icustay_id order by he_time) as next_status,
    LEAD(windowendct, 1) over (partition by SUBJECT_ID, ICUSTAY_ID order by windowendct) as NEXT_ENDCT
  from allevents3
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
        then NEXT_ENDCT else null 
    end he_offset
  from allevents35
)
select h.*, (h.he_offset - h.he_onset) as he_length into he_set 
from hetable1 h where he_onset is not null order by h.subject_id, h.he_onset;

