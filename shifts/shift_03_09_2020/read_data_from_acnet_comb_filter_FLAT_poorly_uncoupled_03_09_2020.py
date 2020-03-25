from datetime import datetime
import fur.path_assistant as path_assistant
from acnet_reader.fur_data_reader import save_acnet_data_for_fur
shift_03_09_2020 = path_assistant.PathAssistant('shift_03_09_2020')
t1 = datetime(2020, 3, 9, 20, 16, 0)
t2 = datetime(2020, 3, 9, 21, 4, 0)
save_acnet_data_for_fur(shift_03_09_2020, t1, t2, "all_acnet_data_comb_filter_FLAT_poorly_uncoupled_03_09_2020.csv")