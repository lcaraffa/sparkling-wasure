## Distributed Surface reconstruction Service


### Algo parameters (PARAMS)
```xml
<env>
  <datasets>
    <austin>
      <plot_lvl>3</plot_lvl>
      <do_process>true</do_process>
      <datatype>files</datatype>
      <dim>2</dim>
      <ndtree_depth>3</ndtree_depth>
      <bbox>0000x10000:0000x10000</bbox>
      <max_ppt>12000</max_ppt>
      <pscale>20</pscale>

      <nb_samples>1</nb_samples>

      <lambda>1</lambda>
      <lambda>10</lambda>

      <mode>surface</mode>
      <algo_seed>1000</algo_seed>
      <StorageLevel>MEMORY_AND_DISK</StorageLevel>
      <!-- Meta parameters -->
      <do_process>true</do_process>
      <do_expand>true</do_expand>
      
    </austin>
  </datasets>
</env>
```