import torch.nn as nn
import torch.nn.functional as F
import torch
import pandas as pd

class Net_59(nn.Module):
    def __init__(self):
        super().__init__()
        self.conv1 = nn.Conv1d(1, 4, 2)
        self.bn1   = nn.BatchNorm1d(4)

        self.conv2 = nn.Conv1d(4, 16, 2)
        self.bn2   = nn.BatchNorm1d(16)

        self.conv3 = nn.Conv1d(16, 32, 2)
        self.bn3   = nn.BatchNorm1d(32)
        
        self.pool = nn.MaxPool1d(2)
        self.fc1 = nn.Linear(32 * 6, 256)
        self.fc2 = nn.Linear(256, 512)
        self.fc3 = nn.Linear(512, 1024)
        self.fc4 = nn.Linear(1024, 59)

    def forward(self, x):
        x = self.pool(F.relu(self.bn1(self.conv1(x))))
        x = self.pool(F.relu(self.bn2(self.conv2(x))))
        x = self.pool(F.relu(self.bn3(self.conv3(x))))
        x = torch.flatten(x, 1) # flatten all dimensions except batch
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        x = F.relu(self.fc3(x))
        # x = F.sigmoid(self.fc4(x))
        x = self.fc4(x)
        return x
    
    
    
def predict_vector(text_id, data_path, df_initial_vector):
    
    model = Net_59()
    model.load_state_dict(torch.load(f"{data_path}/cnn_models/simple_model_trained_atel_gcn_follow_up_60E_20BS_11_09_2024", map_location=torch.device('cpu')))
    model.eval()

    input_tensor  = torch.Tensor(df_initial_vector[text_id])
    in_ = input_tensor.unsqueeze(0).unsqueeze(0)
    out = model(in_).squeeze().detach().numpy()
    out[out<0.01] = 0
    
    dict_out = {}
    dict_out["Legend"] = df_initial_vector["Legend"].values
    dict_out[text_id] = df_initial_vector[text_id].values
    dict_out["Follow-up Vector Prediction"] = out
    return pd.DataFrame(dict_out)
