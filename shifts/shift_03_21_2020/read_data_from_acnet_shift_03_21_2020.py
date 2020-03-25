from datetime import datetime
import fur.path_assistant as path_assistant
from acnet_reader.fur_data_reader import save_acnet_data_for_fur
shift_03_21_2020 = path_assistant.PathAssistant('shift_03_21_2020')
t1 = datetime(2020, 3, 21, 12, 0, 0)
t2 = datetime(2020, 3, 21, 14, 0, 0)
save_acnet_data_for_fur(shift_03_21_2020, t1, t2, "all_acnet_data_for_shift_03_21_2020.csv")