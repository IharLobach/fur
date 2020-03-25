from datetime import datetime
import fur.path_assistant as path_assistant
from acnet_reader.fur_data_reader import save_acnet_data_for_fur
shift_03_16_2020 = path_assistant.PathAssistant('shift_03_16_2020')
t1 = datetime(2020, 3, 16, 13, 4, 0)
t2 = datetime(2020, 3, 16, 13, 7, 0)
save_acnet_data_for_fur(shift_03_16_2020, t1, t2, "all_acnet_data_directly_connected_for_transmission_function_03_16_2020.csv")